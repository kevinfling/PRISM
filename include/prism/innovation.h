/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_INNOVATION_H
#define PRISM_INNOVATION_H
#include <prism.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Borrowed caller-owned state and callback. state_bytes describes the live
 * state region for overlap checks; all reachable state must also be disjoint
 * from outputs. The callback fills count doubles, returns a PRISM status, and
 * must not allocate. Its state and code must outlive every call. To obtain the
 * Gaussian covariance construction, supply independent standard normals;
 * finite-value checks cannot establish their distribution or independence. */
typedef struct {
    void *state;
    size_t state_bytes;
    int (*fill_normal)(void *state, double *out, size_t count, size_t capacity);
} prism_normal_source;
/* count >= 1; capacities count doubles. Descriptor, state, and output must be
 * disjoint. Argument/capacity failures do not invoke the callback. Callback or
 * nonfinite-output failure can consume state and leave partial output. */
int prism_normal_fill(const prism_normal_source *, double *out, size_t count, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
