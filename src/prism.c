/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "internal.h"
#include <float.h>
#include <string.h>
struct prism_plan {
    size_t n, count, bytes, component_values, weight_values;
    prism_allocator allocator;
    const prism_backend *backend;
    double *weights;
};
struct prism_workspace {
    const prism_plan *plan;
    double *spectrum, *filtered;
    void *backend_state;
};
static void *default_allocate(void *user, size_t bytes) {
    (void)user;
    return malloc(bytes);
}
static void default_deallocate(void *user, void *p) {
    (void)user;
    free(p);
}
static int finite_array(const double *x, size_t n) {
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(x[i]))
            return 0;
    return 1;
}
int prism_plan_create(const prism_config *c, const prism_backend *b,
                      const prism_allocator *allocator, prism_plan **out) {
    size_t bytes, values, weight_values, weight_bytes;
    prism_allocator a = {NULL, default_allocate, default_deallocate};
    if (!out)
        return PRISM_EINVAL;
    *out = NULL;
    if (!c || c->n < 2 || !c->component_count || !isfinite(c->sample_spacing) ||
        c->sample_spacing <= 0)
        return PRISM_EINVAL;
    if (c->component_count == SIZE_MAX || !prism_size_mul(c->n, sizeof(double), &bytes) ||
        !prism_size_mul(c->n, c->component_count, &values) ||
        !prism_size_mul(c->n, c->component_count + 1, &weight_values) ||
        !prism_size_mul(weight_values, sizeof(double), &weight_bytes))
        return PRISM_EOVERFLOW;
    long double length = (long double)c->n * (long double)c->sample_spacing;
    if (length > (long double)DBL_MAX)
        return PRISM_ENUMERIC;
    if (allocator) {
        if (!allocator->allocate || !allocator->deallocate)
            return PRISM_EINVAL;
        a = *allocator;
    }
    if (!b)
        b = prism_backend_fftw();
    if (!b)
        b = prism_backend_reference();
    if (!b->name || !b->create || !b->forward || !b->inverse || !b->destroy)
        return PRISM_EINVAL;
    if (c->mode == PRISM_BANDS_FIXED) {
        if (!c->band_edges || c->band_edges[0] != c->n)
            return PRISM_EINVAL;
        for (size_t j = 0; j < c->component_count; ++j)
            if (!c->band_edges[j + 1] || c->band_edges[j + 1] >= c->band_edges[j])
                return PRISM_EINVAL;
    } else if (c->mode == PRISM_BANDS_HEAT_TELESCOPE) {
        if (!c->heat_times)
            return PRISM_EINVAL;
        double prev = 0;
        for (size_t j = 0; j < c->component_count; ++j) {
            if (!isfinite(c->heat_times[j]) || c->heat_times[j] <= prev)
                return PRISM_EINVAL;
            prev = c->heat_times[j];
        }
    } else
        return PRISM_EINVAL;
    prism_plan *p = a.allocate(a.user, sizeof(*p));
    if (!p)
        return PRISM_ENOMEM;
    p->allocator = a;
    p->n = c->n;
    p->count = c->component_count;
    p->bytes = bytes;
    p->backend = b;
    p->component_values = values;
    p->weight_values = weight_values;
    p->weights = a.allocate(a.user, weight_bytes);
    if (!p->weights) {
        a.deallocate(a.user, p);
        return PRISM_ENOMEM;
    }
    if (c->mode == PRISM_BANDS_FIXED) {
        for (size_t j = 0; j <= p->count; ++j)
            for (size_t k = 0; k < p->n; ++k)
                p->weights[j * p->n + k] =
                    (j == p->count ? k < c->band_edges[j]
                                   : k >= c->band_edges[j + 1] && k < c->band_edges[j])
                        ? 1.0
                        : 0.0;
    } else {
        for (size_t k = 0; k < p->n; ++k) {
            long double omega = PRISM_PI * (long double)k / length;
            long double lambda = omega * omega, prev = 1.0L, prev_time = 0;
            for (size_t j = 0; j < p->count; ++j) {
                long double t = (long double)c->heat_times[j];
                long double next = k ? expl(-t * lambda) : 1.0L;
                long double difference = k ? prev * (-expm1l(-(t - prev_time) * lambda)) : 0.0L;
                p->weights[j * p->n + k] = (double)difference;
                prev = next;
                prev_time = t;
            }
            p->weights[p->count * p->n + k] = (double)prev;
        }
    }
    if (!finite_array(p->weights, weight_values)) {
        prism_plan_destroy(p);
        return PRISM_ENUMERIC;
    }
    *out = p;
    return PRISM_OK;
}
void prism_plan_destroy(prism_plan *p) {
    if (!p)
        return;
    prism_allocator a = p->allocator;
    a.deallocate(a.user, p->weights);
    a.deallocate(a.user, p);
}
int prism_workspace_create(const prism_plan *p, prism_workspace **out) {
    if (!out)
        return PRISM_EINVAL;
    *out = NULL;
    if (!p)
        return PRISM_EINVAL;
    prism_allocator a = p->allocator;
    prism_workspace *w = a.allocate(a.user, sizeof(*w));
    if (!w)
        return PRISM_ENOMEM;
    memset(w, 0, sizeof(*w));
    w->plan = p;
    w->spectrum = a.allocate(a.user, p->bytes);
    w->filtered = a.allocate(a.user, p->bytes);
    if (!w->spectrum || !w->filtered) {
        prism_workspace_destroy(w);
        return PRISM_ENOMEM;
    }
    int status = p->backend->create(p->n, &a, &w->backend_state);
    if (status != PRISM_OK) {
        prism_workspace_destroy(w);
        return status;
    }
    if (!w->backend_state) {
        prism_workspace_destroy(w);
        return PRISM_EBACKEND;
    }
    *out = w;
    return PRISM_OK;
}
void prism_workspace_destroy(prism_workspace *w) {
    if (!w)
        return;
    prism_allocator a = w->plan->allocator;
    if (w->backend_state)
        w->plan->backend->destroy(w->backend_state);
    if (w->spectrum)
        a.deallocate(a.user, w->spectrum);
    if (w->filtered)
        a.deallocate(a.user, w->filtered);
    a.deallocate(a.user, w);
}
int prism_plan_weights(const prism_plan *p, double *out, size_t capacity) {
    if (!p || !out)
        return PRISM_EINVAL;
    if (capacity < p->weight_values)
        return PRISM_ECAPACITY;
    memcpy(out, p->weights, p->weight_values * sizeof(double));
    return PRISM_OK;
}
int prism_decompose(const prism_plan *p, prism_workspace *w, const double *x, prism_result *out) {
    if (!out)
        return PRISM_EINVAL;
    out->count = 0;
    out->n = 0;
    out->reconstruction_relative_error = 0;
    out->reconstruction_max_abs_error = 0;
    if (!p || !w || w->plan != p || !x || !out->components || !out->residual)
        return PRISM_EINVAL;
    if (out->component_capacity < p->component_values || out->residual_capacity < p->n)
        return PRISM_ECAPACITY;
    size_t component_bytes = p->component_values * sizeof(double);
    if (prism_overlap(x, p->bytes, out->components, component_bytes) ||
        prism_overlap(x, p->bytes, out->residual, p->bytes) ||
        prism_overlap(out->components, component_bytes, out->residual, p->bytes))
        return PRISM_EINVAL;
    if (!finite_array(x, p->n))
        return PRISM_ENUMERIC;
    int status = p->backend->forward(w->backend_state, x, w->spectrum);
    if (status != PRISM_OK)
        return status;
    if (!finite_array(w->spectrum, p->n))
        return PRISM_ENUMERIC;
    for (size_t j = 0; j <= p->count; ++j) {
        for (size_t k = 0; k < p->n; ++k)
            w->filtered[k] = w->spectrum[k] * p->weights[j * p->n + k];
        double *y = j == p->count ? out->residual : out->components + j * p->n;
        status = p->backend->inverse(w->backend_state, w->filtered, y);
        if (status != PRISM_OK)
            return status;
        if (!finite_array(y, p->n))
            return PRISM_ENUMERIC;
    }
    long double error = 0, norm = 0, max_error = 0;
    for (size_t i = 0; i < p->n; ++i) {
        long double reconstructed = (long double)out->residual[i];
        for (size_t j = 0; j < p->count; ++j)
            reconstructed += (long double)out->components[j * p->n + i];
        long double difference = fabsl((long double)x[i] - reconstructed);
        error = hypotl(error, difference);
        norm = hypotl(norm, (long double)x[i]);
        if (difference > max_error)
            max_error = difference;
    }
    long double relative = norm > 0 ? error / norm : error;
    if (!isfinite(error) || !isfinite(norm) || !isfinite(max_error) || !isfinite(relative) ||
        max_error > (long double)DBL_MAX || relative > (long double)DBL_MAX)
        return PRISM_ENUMERIC;
    out->reconstruction_relative_error = (double)relative;
    out->reconstruction_max_abs_error = (double)max_error;
    out->count = p->count;
    out->n = p->n;
    return PRISM_OK;
}
int prism_reconstruct(const prism_result *r, double *out, size_t capacity) {
    size_t bytes, values, component_bytes;
    if (!r || !out || !r->n || !r->count || !r->components || !r->residual)
        return PRISM_EINVAL;
    if (!prism_size_mul(r->n, sizeof(double), &bytes) || !prism_size_mul(r->n, r->count, &values) ||
        !prism_size_mul(values, sizeof(double), &component_bytes))
        return PRISM_EOVERFLOW;
    if (capacity < r->n || r->component_capacity < values || r->residual_capacity < r->n)
        return PRISM_ECAPACITY;
    if (prism_overlap(out, bytes, r->components, component_bytes) ||
        prism_overlap(out, bytes, r->residual, bytes) ||
        prism_overlap(r->components, component_bytes, r->residual, bytes))
        return PRISM_EINVAL;
    if (!finite_array(r->components, values) || !finite_array(r->residual, r->n))
        return PRISM_ENUMERIC;
    for (size_t i = 0; i < r->n; ++i) {
        long double sum = (long double)r->residual[i];
        for (size_t j = 0; j < r->count; ++j)
            sum += (long double)r->components[j * r->n + i];
        if (!isfinite(sum) || fabsl(sum) > (long double)DBL_MAX)
            return PRISM_ENUMERIC;
        out[i] = (double)sum;
    }
    return PRISM_OK;
}
const char *prism_status_string(int status) {
    switch (status) {
    case PRISM_OK:
        return "success";
    case PRISM_EINVAL:
        return "invalid argument or overlapping buffers";
    case PRISM_ECAPACITY:
        return "insufficient output capacity";
    case PRISM_ENOMEM:
        return "allocation failed";
    case PRISM_EBACKEND:
        return "transform backend failed";
    case PRISM_ENUMERIC:
        return "nonfinite data or numerical range exceeded";
    case PRISM_EOVERFLOW:
        return "size exceeds supported range";
    case PRISM_ECOVARIANCE:
        return "covariance factorization or circulant embedding rejected";
    default:
        return "unknown status";
    }
}
