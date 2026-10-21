/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "internal.h"
typedef struct {
    size_t n;
    prism_allocator allocator;
} reference_state;
static int reference_create(size_t n, const prism_allocator *a, void **out) {
    reference_state *s = a->allocate(a->user, sizeof(*s));
    if (!s)
        return PRISM_ENOMEM;
    s->n = n;
    s->allocator = *a;
    *out = s;
    return PRISM_OK;
}
static int reference_forward(void *state, const double *x, double *y) {
    const reference_state *s = state;
    for (size_t k = 0; k < s->n; ++k) {
        double sum = 0, correction = 0;
        double angle = (double)PRISM_PI * (double)k / (double)s->n;
        double dc = cos(angle), ds = sin(angle), c = 0, v = 0;
        for (size_t i = 0; i < s->n; ++i) {
            /* Re-anchor rotation to limit drift on long vectors. */
            if (i % 128 == 0) {
                c = cos(((double)i + 0.5) * angle);
                v = sin(((double)i + 0.5) * angle);
            }
            double term = x[i] * c - correction;
            double added = sum + term;
            correction = (added - sum) - term;
            sum = added;
            double next = c * dc - v * ds;
            v = v * dc + c * ds;
            c = next;
        }
        y[k] = (double)(sum * sqrt((k ? 2.0 : 1.0) / (double)s->n));
    }
    return PRISM_OK;
}
static int reference_inverse(void *state, const double *x, double *y) {
    const reference_state *s = state;
    double scale = sqrt(2.0 / (double)s->n);
    for (size_t i = 0; i < s->n; ++i) {
        double sum = x[0] / sqrt((double)s->n), correction = 0;
        double angle = (double)PRISM_PI * ((double)i + 0.5) / (double)s->n;
        double dc = cos(angle), ds = sin(angle), c = dc, v = ds;
        for (size_t k = 1; k < s->n; ++k) {
            if (k % 128 == 0) {
                c = cos((double)k * angle);
                v = sin((double)k * angle);
            }
            double term = scale * x[k] * c - correction;
            double added = sum + term;
            correction = (added - sum) - term;
            sum = added;
            double next = c * dc - v * ds;
            v = v * dc + c * ds;
            c = next;
        }
        y[i] = (double)sum;
    }
    return PRISM_OK;
}
static void reference_destroy(void *state) {
    reference_state *s = state;
    prism_allocator a = s->allocator;
    a.deallocate(a.user, s);
}
const prism_backend *prism_backend_reference(void) {
    static const prism_backend b = {"reference", reference_create, reference_forward,
                                    reference_inverse, reference_destroy};
    return &b;
}
