/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "internal.h"
#include <float.h>
#include <limits.h>
#include <prism/noise.h>
#include <string.h>
#ifdef PRISM_HAVE_FFTW
#include <fftw3.h>
#endif

struct prism_noise_plan {
    prism_allocator allocator;
    prism_noise_method method;
    size_t n, m;
    double *lags, *factor;
};
struct prism_noise_workspace {
    const prism_noise_plan *owner;
    double *white, *candidate, *spectrum;
#ifdef PRISM_HAVE_FFTW
    fftw_plan forward, inverse;
#endif
};
static void *default_allocate(void *user, size_t bytes) {
    (void)user;
    return malloc(bytes);
}
static void default_deallocate(void *user, void *ptr) {
    (void)user;
    free(ptr);
}
static void *allocate(const prism_noise_plan *p, size_t count) {
    return p->allocator.allocate(p->allocator.user, count * sizeof(double));
}
static void release(const prism_noise_plan *p, void *ptr) {
    if (ptr)
        p->allocator.deallocate(p->allocator.user, ptr);
}
int prism_noise_method_available(prism_noise_method method) {
    if (method == PRISM_NOISE_CHOLESKY || method == PRISM_NOISE_CIRCULANT_REFERENCE)
        return 1;
#ifdef PRISM_HAVE_FFTW
    if (method == PRISM_NOISE_CIRCULANT_FFTW)
        return 1;
#endif
    return 0;
}

int prism_fgn_covariance(size_t count, double hurst, double variance, double *out,
                         size_t capacity) {
    size_t bytes;
    if (!count || !out || !isfinite(hurst) || hurst <= 0 || hurst >= 1 || !isfinite(variance) ||
        variance < 0)
        return PRISM_EINVAL;
    if (!prism_size_mul(count, sizeof(double), &bytes))
        return PRISM_EOVERFLOW;
    if (capacity < count)
        return PRISM_ECAPACITY;
    out[0] = variance;
    long double p = 2 * (long double)hurst;
    for (size_t k = 1; k < count; ++k) {
        long double gamma;
        if (hurst == 0.5 || variance == 0) {
            gamma = 0;
        } else if (k == 1) {
            gamma = expm1l((p - 1) * logl(2));
        } else {
            /* Even binomial expansion of the second difference. Avoids
             * subtracting nearly equal powers at large lags or H near 1/2. */
            long double lag = (long double)k;
            long double inv_square = 1 / (lag * lag);
            long double term = 1, sum = 1;
            for (unsigned j = 1; j < 128; ++j) {
                long double a = 2 * (long double)j;
                term *= (p - a) * (p - a - 1) / ((a + 1) * (a + 2)) * inv_square;
                sum += term;
                if (fabsl(term) <= LDBL_EPSILON * fabsl(sum))
                    break;
            }
            gamma = 0.5L * p * (p - 1) * expl((p - 2) * logl(lag)) * sum;
        }
        out[k] = (double)((long double)variance * gamma);
    }
    return PRISM_OK;
}

/* Real DFT, unnormalized forward and inverse divided by m. These deliberately
 * direct O(m^2) transforms are the independent, dependency-free reference. */
static void forward_reference(prism_noise_workspace *w) {
    size_t m = w->owner->m;
    for (size_t k = 0; k <= m / 2; ++k) {
        long double re = 0, im = 0;
        for (size_t j = 0; j < m; ++j) {
            long double angle = 2 * PRISM_PI * (long double)k * (long double)j / (long double)m;
            re += (long double)w->white[j] * cosl(angle);
            im -= (long double)w->white[j] * sinl(angle);
        }
        w->spectrum[2 * k] = (double)re;
        w->spectrum[2 * k + 1] = (double)im;
    }
}
static void inverse_reference(prism_noise_workspace *w) {
    size_t m = w->owner->m;
    for (size_t j = 0; j < m; ++j) {
        long double sum = w->spectrum[0];
        for (size_t k = 1; k <= m / 2; ++k) {
            long double angle = 2 * PRISM_PI * (long double)k * (long double)j / (long double)m;
            long double contribution = (long double)w->spectrum[2 * k] * cosl(angle) -
                                       (long double)w->spectrum[2 * k + 1] * sinl(angle);
            sum += (m % 2 == 0 && k == m / 2) ? contribution : 2 * contribution;
        }
        w->candidate[j] = (double)(sum / (long double)m);
    }
}
static void forward(prism_noise_workspace *w) {
#ifdef PRISM_HAVE_FFTW
    if (w->owner->method == PRISM_NOISE_CIRCULANT_FFTW) {
        fftw_execute(w->forward);
        return;
    }
#endif
    forward_reference(w);
}
static void inverse(prism_noise_workspace *w) {
#ifdef PRISM_HAVE_FFTW
    if (w->owner->method == PRISM_NOISE_CIRCULANT_FFTW) {
        fftw_execute(w->inverse);
        for (size_t i = 0; i < w->owner->m; ++i)
            w->candidate[i] /= (double)w->owner->m;
        return;
    }
#endif
    inverse_reference(w);
}

void prism_noise_workspace_destroy(prism_noise_workspace *w) {
    if (!w)
        return;
    const prism_noise_plan *p = w->owner;
#ifdef PRISM_HAVE_FFTW
    if (w->forward)
        fftw_destroy_plan(w->forward);
    if (w->inverse)
        fftw_destroy_plan(w->inverse);
#endif
    release(p, w->spectrum);
    release(p, w->candidate);
    release(p, w->white);
    release(p, w);
}
int prism_noise_workspace_create(const prism_noise_plan *p, prism_noise_workspace **out) {
    if (!out)
        return PRISM_EINVAL;
    *out = NULL;
    if (!p)
        return PRISM_EINVAL;
    prism_noise_workspace *w = p->allocator.allocate(p->allocator.user, sizeof(*w));
    if (!w)
        return PRISM_ENOMEM;
    memset(w, 0, sizeof(*w));
    w->owner = p;
    w->white = allocate(p, p->m);
    w->candidate = allocate(p, p->m);
    if (p->method != PRISM_NOISE_CHOLESKY)
        w->spectrum = allocate(p, 2 * (p->m / 2 + 1));
    if (!w->white || !w->candidate || (p->method != PRISM_NOISE_CHOLESKY && !w->spectrum)) {
        prism_noise_workspace_destroy(w);
        return PRISM_ENOMEM;
    }
#ifdef PRISM_HAVE_FFTW
    if (p->method == PRISM_NOISE_CIRCULANT_FFTW) {
        w->forward = fftw_plan_dft_r2c_1d((int)p->m, w->white, (fftw_complex *)(void *)w->spectrum,
                                          FFTW_ESTIMATE | FFTW_UNALIGNED);
        w->inverse = fftw_plan_dft_c2r_1d((int)p->m, (fftw_complex *)(void *)w->spectrum,
                                          w->candidate, FFTW_ESTIMATE | FFTW_UNALIGNED);
        if (!w->forward || !w->inverse) {
            prism_noise_workspace_destroy(w);
            return PRISM_EBACKEND;
        }
    }
#endif
    *out = w;
    return PRISM_OK;
}
void prism_noise_plan_destroy(prism_noise_plan *p) {
    if (!p)
        return;
    release(p, p->factor);
    release(p, p->lags);
    prism_allocator a = p->allocator;
    a.deallocate(a.user, p);
}
static int cholesky(prism_noise_plan *p, prism_noise_diagnostics *d) {
    memset(p->factor, 0, p->n * p->n * sizeof(double));
    d->minimum_cholesky_pivot = p->lags[0];
    for (size_t i = 0; i < p->n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            long double sum = p->lags[i - j];
            for (size_t k = 0; k < j; ++k)
                sum -= (long double)p->factor[i * p->n + k] * p->factor[j * p->n + k];
            if (!isfinite(sum))
                return PRISM_ENUMERIC;
            if (i == j) {
                if (sum < d->minimum_cholesky_pivot)
                    d->minimum_cholesky_pivot = (double)sum;
                if (sum < 0)
                    return PRISM_ECOVARIANCE;
                p->factor[i * p->n + j] = (double)sqrtl(sum);
            } else if (p->factor[j * p->n + j] == 0) {
                if (sum != 0)
                    return PRISM_ECOVARIANCE;
            } else {
                p->factor[i * p->n + j] = (double)(sum / p->factor[j * p->n + j]);
                if (!isfinite(p->factor[i * p->n + j]))
                    return PRISM_ENUMERIC;
            }
        }
    }
    return PRISM_OK;
}
static int embedding(prism_noise_plan *p, double tolerance, prism_noise_diagnostics *d) {
    prism_noise_workspace *w = NULL;
    int status = prism_noise_workspace_create(p, &w);
    if (status != PRISM_OK)
        return status;
    for (size_t j = 0; j < p->m; ++j)
        w->white[j] = p->lags[j < p->n ? j : p->m - j];
    forward(w);
    d->embedding_size = p->m;
    d->minimum_embedding_eigenvalue = w->spectrum[0];
    d->maximum_embedding_eigenvalue = w->spectrum[0];
    double max_absolute = 0;
    for (size_t k = 0; k <= p->m / 2; ++k) {
        double re = w->spectrum[2 * k], im = w->spectrum[2 * k + 1];
        if (!isfinite(re) || !isfinite(im)) {
            status = PRISM_ENUMERIC;
            goto done;
        }
        d->minimum_embedding_eigenvalue = fmin(d->minimum_embedding_eigenvalue, re);
        d->maximum_embedding_eigenvalue = fmax(d->maximum_embedding_eigenvalue, re);
        d->maximum_embedding_imaginary_part = fmax(d->maximum_embedding_imaginary_part, fabs(im));
        max_absolute = fmax(max_absolute, fabs(re));
    }
    for (size_t k = 0; k <= p->m / 2; ++k) {
        double eigenvalue = w->spectrum[2 * k];
        if (eigenvalue < 0) {
            if (-eigenvalue > tolerance * max_absolute) {
                status = PRISM_ECOVARIANCE;
                goto done;
            }
            d->rounded_eigenvalue_count += (k == 0 || (p->m % 2 == 0 && k == p->m / 2)) ? 1u : 2u;
            d->maximum_eigenvalue_adjustment = fmax(d->maximum_eigenvalue_adjustment, -eigenvalue);
            eigenvalue = 0;
        }
        p->factor[k] = sqrt(eigenvalue);
    }
done:
    prism_noise_workspace_destroy(w);
    return status;
}
int prism_noise_plan_create(const prism_noise_config *config, const double *lags, size_t capacity,
                            const prism_allocator *allocator, prism_noise_plan **out,
                            prism_noise_diagnostics *d) {
    if (out)
        *out = NULL;
    if (d)
        memset(d, 0, sizeof(*d));
    if (!out || !d || !config || !lags || !config->n ||
        !isfinite(config->embedding_relative_tolerance) ||
        config->embedding_relative_tolerance < 0 || config->embedding_relative_tolerance >= 1)
        return PRISM_EINVAL;
    d->method = config->method;
    d->n = config->n;
    if (!prism_noise_method_available(config->method))
        return config->method == PRISM_NOISE_CIRCULANT_FFTW ? PRISM_EBACKEND : PRISM_EINVAL;
    if (config->method == PRISM_NOISE_CHOLESKY && config->embedding_relative_tolerance != 0)
        return PRISM_EINVAL;
    prism_allocator a = {NULL, default_allocate, default_deallocate};
    if (allocator) {
        if (!allocator->allocate || !allocator->deallocate)
            return PRISM_EINVAL;
        a = *allocator;
    }
    size_t n = config->n, m = n, factor_count, bytes;
    if (config->method == PRISM_NOISE_CHOLESKY) {
        if (!prism_size_mul(n, n, &factor_count))
            return PRISM_EOVERFLOW;
    } else {
        if (n > 1 && !prism_size_mul(n - 1, 2, &m))
            return PRISM_EOVERFLOW;
        if (m / 2 + 1 > SIZE_MAX / 2)
            return PRISM_EOVERFLOW;
        factor_count = m / 2 + 1;
        if (!prism_size_mul(2 * factor_count, sizeof(double), &bytes))
            return PRISM_EOVERFLOW;
        if (config->method == PRISM_NOISE_CIRCULANT_FFTW && m > INT_MAX)
            return PRISM_EOVERFLOW;
    }
    if (!prism_size_mul(factor_count, sizeof(double), &bytes) ||
        !prism_size_mul(m, sizeof(double), &bytes) || !prism_size_mul(n, sizeof(double), &bytes))
        return PRISM_EOVERFLOW;
    if (capacity < n)
        return PRISM_ECAPACITY;
    for (size_t i = 0; i < n; ++i)
        if (!isfinite(lags[i]))
            return PRISM_EINVAL;
    if (lags[0] < 0)
        return PRISM_ECOVARIANCE;
    for (size_t i = 1; i < n; ++i)
        if (fabs(lags[i]) > lags[0])
            return PRISM_ECOVARIANCE;
    prism_noise_plan *p = a.allocate(a.user, sizeof(*p));
    if (!p)
        return PRISM_ENOMEM;
    *p = (prism_noise_plan){a, config->method, n, m, NULL, NULL};
    p->lags = allocate(p, n);
    p->factor = allocate(p, factor_count);
    if (!p->lags || !p->factor) {
        prism_noise_plan_destroy(p);
        return PRISM_ENOMEM;
    }
    memcpy(p->lags, lags, n * sizeof(double));
    int status = config->method == PRISM_NOISE_CHOLESKY
                     ? cholesky(p, d)
                     : embedding(p, config->embedding_relative_tolerance, d);
    if (status != PRISM_OK) {
        prism_noise_plan_destroy(p);
        return status;
    }
    *out = p;
    return PRISM_OK;
}
static int apply_covariance(const prism_noise_plan *p, prism_noise_workspace *w, double *out) {
    if (p->method == PRISM_NOISE_CHOLESKY) {
        for (size_t i = 0; i < p->n; ++i) {
            long double sum = 0;
            for (size_t j = 0; j <= i; ++j)
                sum += (long double)p->factor[i * p->n + j] * w->white[j];
            w->candidate[i] = (double)sum;
        }
    } else {
        forward(w);
        for (size_t k = 0; k <= p->m / 2; ++k) {
            w->spectrum[2 * k] *= p->factor[k];
            w->spectrum[2 * k + 1] *= p->factor[k];
        }
        inverse(w);
    }
    for (size_t i = 0; i < p->n; ++i)
        if (!isfinite(w->candidate[i]))
            return PRISM_ENUMERIC;
    memcpy(out, w->candidate, p->n * sizeof(double));
    return PRISM_OK;
}

size_t prism_noise_innovation_count(const prism_noise_plan *p) { return p ? p->m : 0; }
int prism_noise_transform(const prism_noise_plan *p, prism_noise_workspace *w,
                          const double *innovations, size_t innovation_capacity, double *out,
                          size_t capacity) {
    if (!p || !w || w->owner != p || !innovations || !out)
        return PRISM_EINVAL;
    if (capacity < p->n || innovation_capacity < p->m)
        return PRISM_ECAPACITY;
    if (prism_overlap(innovations, p->m * sizeof(double), out, p->n * sizeof(double)))
        return PRISM_EINVAL;
    for (size_t i = 0; i < p->m; ++i)
        if (!isfinite(innovations[i]))
            return PRISM_ENUMERIC;
    memcpy(w->white, innovations, p->m * sizeof(double));
    return apply_covariance(p, w, out);
}
int prism_noise_sample_from(const prism_noise_plan *p, prism_noise_workspace *w,
                            const prism_normal_source *source, double *out, size_t capacity) {
    if (!p || !w || w->owner != p || !source || !out)
        return PRISM_EINVAL;
    if (capacity < p->n)
        return PRISM_ECAPACITY;
    if (prism_overlap(source, sizeof(*source), out, p->n * sizeof(double)) ||
        (source->state &&
         prism_overlap(source->state, source->state_bytes, out, p->n * sizeof(double))))
        return PRISM_EINVAL;
    int status = prism_normal_fill(source, w->white, p->m, p->m);
    if (status != PRISM_OK)
        return status;
    return apply_covariance(p, w, out);
}
int prism_noise_sample(const prism_noise_plan *p, prism_noise_workspace *w, prism_rng *rng,
                       double *out, size_t capacity) {
    prism_normal_source source = prism_rng_source(rng);
    return prism_noise_sample_from(p, w, &source, out, capacity);
}
int prism_noise_target_covariance(const prism_noise_plan *p, double *out, size_t capacity) {
    if (!p || !out)
        return PRISM_EINVAL;
    if (capacity < p->n)
        return PRISM_ECAPACITY;
    memcpy(out, p->lags, p->n * sizeof(double));
    return PRISM_OK;
}
