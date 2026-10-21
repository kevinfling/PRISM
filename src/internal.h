/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_INTERNAL_H
#define PRISM_INTERNAL_H
#include <math.h>
#include <prism.h>
#include <stdint.h>
#include <stdlib.h>
#define PRISM_PI 3.141592653589793238462643383279502884L
static inline int prism_size_mul(size_t a, size_t b, size_t *out) {
    if (b && a > SIZE_MAX / b)
        return 0;
    *out = a * b;
    return 1;
}
static inline int prism_overlap(const void *a, size_t an, const void *b, size_t bn) {
    uintptr_t x = (uintptr_t)a, y = (uintptr_t)b;
    if (an > UINTPTR_MAX - x || bn > UINTPTR_MAX - y)
        return 1;
    return x <= y ? y - x < an : x - y < bn;
}
#endif
