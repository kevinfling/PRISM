/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "internal.h"
#include <prism/innovation.h>

int prism_normal_fill(const prism_normal_source *source, double *out, size_t count,
                      size_t capacity) {
    size_t bytes;
    if (!source || !source->state || !source->state_bytes || !source->fill_normal || !out || !count)
        return PRISM_EINVAL;
    if (!prism_size_mul(count, sizeof(double), &bytes))
        return PRISM_EOVERFLOW;
    if (capacity < count)
        return PRISM_ECAPACITY;
    if (prism_overlap(source, sizeof(*source), out, bytes) ||
        prism_overlap(source->state, source->state_bytes, out, bytes) ||
        prism_overlap(source, sizeof(*source), source->state, source->state_bytes))
        return PRISM_EINVAL;
    int status = source->fill_normal(source->state, out, count, capacity);
    if (status != PRISM_OK)
        return status;
    for (size_t i = 0; i < count; ++i)
        if (!isfinite(out[i]))
            return PRISM_ENUMERIC;
    return PRISM_OK;
}
