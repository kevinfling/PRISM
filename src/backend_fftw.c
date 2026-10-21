/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "internal.h"
#ifdef PRISM_HAVE_FFTW
#include <fftw3.h>
#include <limits.h>
#include <string.h>
typedef struct {
    size_t n;
    prism_allocator allocator;
    double *buffer;
    fftw_plan forward, inverse;
} fftw_state;
static void backend_destroy(void *state) {
    fftw_state *s = state;
    prism_allocator a = s->allocator;
    if (s->forward)
        fftw_destroy_plan(s->forward);
    if (s->inverse)
        fftw_destroy_plan(s->inverse);
    if (s->buffer)
        a.deallocate(a.user, s->buffer);
    a.deallocate(a.user, s);
}
static int backend_create(size_t n, const prism_allocator *a, void **out) {
    size_t bytes;
    if (n > (size_t)INT_MAX || !prism_size_mul(n, sizeof(double), &bytes))
        return PRISM_EOVERFLOW;
    fftw_state *s = a->allocate(a->user, sizeof(*s));
    if (!s)
        return PRISM_ENOMEM;
    memset(s, 0, sizeof(*s));
    s->n = n;
    s->allocator = *a;
    s->buffer = a->allocate(a->user, bytes);
    if (!s->buffer) {
        backend_destroy(s);
        return PRISM_ENOMEM;
    }
    /* UNALIGNED supports ordinary and caller-provided allocators. */
    s->forward = fftw_plan_r2r_1d((int)n, s->buffer, s->buffer, FFTW_REDFT10,
                                  FFTW_ESTIMATE | FFTW_UNALIGNED);
    s->inverse = fftw_plan_r2r_1d((int)n, s->buffer, s->buffer, FFTW_REDFT01,
                                  FFTW_ESTIMATE | FFTW_UNALIGNED);
    if (!s->forward || !s->inverse) {
        backend_destroy(s);
        return PRISM_EBACKEND;
    }
    *out = s;
    return PRISM_OK;
}
static int backend_forward(void *state, const double *x, double *y) {
    fftw_state *s = state;
    memcpy(s->buffer, x, s->n * sizeof(double));
    fftw_execute(s->forward);
    double scale = sqrt(0.5 / (double)s->n);
    y[0] = s->buffer[0] * 0.5 / sqrt((double)s->n);
    for (size_t k = 1; k < s->n; ++k)
        y[k] = s->buffer[k] * scale;
    return PRISM_OK;
}
static int backend_inverse(void *state, const double *x, double *y) {
    fftw_state *s = state;
    s->buffer[0] = x[0] / sqrt((double)s->n);
    double scale = sqrt(0.5 / (double)s->n);
    for (size_t k = 1; k < s->n; ++k)
        s->buffer[k] = x[k] * scale;
    fftw_execute(s->inverse);
    memcpy(y, s->buffer, s->n * sizeof(double));
    return PRISM_OK;
}
const prism_backend *prism_backend_fftw(void) {
    static const prism_backend b = {"fftw3", backend_create, backend_forward, backend_inverse,
                                    backend_destroy};
    return &b;
}
const char *prism_fftw_version(void) { return fftw_version; }
#else
const prism_backend *prism_backend_fftw(void) { return NULL; }
const char *prism_fftw_version(void) { return "disabled"; }
#endif
