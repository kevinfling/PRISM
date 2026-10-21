/* SPDX-License-Identifier: Apache-2.0
 * PCG Random Number Generation for C.
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 * Copyright (c) 2026 Kevin Fling
 *
 * PCG transition, output permutation, and seeding adapted from pcg-c-basic.
 * Modified for PRISM: checked caller-owned API, explicit sequencing, 53-bit
 * uniforms, and cached Box-Muller normals. See docs/STOCHASTIC.md.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *     https://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. A copy is in LICENSES/Apache-2.0.txt. */
#include "internal.h"
#include <prism/rng.h>

static uint32_t step(prism_rng *rng) {
    uint64_t previous = rng->state;
    rng->state = previous * UINT64_C(6364136223846793005) + rng->increment;
    uint32_t mixed = (uint32_t)(((previous >> 18u) ^ previous) >> 27u);
    uint32_t rotation = (uint32_t)(previous >> 59u);
    return (mixed >> rotation) | (mixed << ((32u - rotation) & 31u));
}

static int valid(const prism_rng *rng) {
    return rng && (rng->increment & UINT64_C(1)) && rng->has_normal_spare <= 1u &&
           (!rng->has_normal_spare || isfinite(rng->normal_spare));
}

static double uniform(prism_rng *rng) {
    uint64_t high = step(rng);
    uint64_t low = step(rng);
    uint64_t bits = (high << 32u) | low;
    return (double)(bits >> 11u) * 0x1p-53;
}

int prism_rng_seed(prism_rng *rng, uint64_t seed, uint64_t stream) {
    if (!rng || stream > (UINT64_MAX >> 1u))
        return PRISM_EINVAL;
    rng->state = 0;
    rng->increment = (stream << 1u) | UINT64_C(1);
    rng->normal_spare = 0;
    rng->has_normal_spare = 0;
    (void)step(rng);
    rng->state += seed;
    (void)step(rng);
    return PRISM_OK;
}

int prism_rng_next_u32(prism_rng *rng, uint32_t *out) {
    if (!valid(rng) || !out || prism_overlap(rng, sizeof(*rng), out, sizeof(*out)))
        return PRISM_EINVAL;
    *out = step(rng);
    return PRISM_OK;
}

static int validate_buffer(prism_rng *rng, double *out, size_t count, size_t capacity) {
    size_t bytes;
    if (!valid(rng) || !out || !count)
        return PRISM_EINVAL;
    if (!prism_size_mul(count, sizeof(double), &bytes))
        return PRISM_EOVERFLOW;
    if (capacity < count)
        return PRISM_ECAPACITY;
    if (prism_overlap(rng, sizeof(*rng), out, bytes))
        return PRISM_EINVAL;
    return PRISM_OK;
}

int prism_rng_uniform(prism_rng *rng, double *out, size_t count, size_t capacity) {
    int status = validate_buffer(rng, out, count, capacity);
    if (status != PRISM_OK)
        return status;
    for (size_t i = 0; i < count; ++i)
        out[i] = uniform(rng);
    return PRISM_OK;
}

int prism_rng_normal(prism_rng *rng, double *out, size_t count, size_t capacity) {
    int status = validate_buffer(rng, out, count, capacity);
    if (status != PRISM_OK)
        return status;
    for (size_t i = 0; i < count; ++i) {
        if (rng->has_normal_spare) {
            out[i] = rng->normal_spare;
            rng->has_normal_spare = 0;
        } else {
            /* 1-U is in (0,1], including when U is zero. */
            double radius = sqrt(-2 * log(1 - uniform(rng)));
            double angle = 2 * (double)PRISM_PI * uniform(rng);
            out[i] = radius * cos(angle);
            rng->normal_spare = radius * sin(angle);
            rng->has_normal_spare = 1;
        }
    }
    return PRISM_OK;
}

static int fill_normal(void *state, double *out, size_t count, size_t capacity) {
    return prism_rng_normal(state, out, count, capacity);
}
prism_normal_source prism_rng_source(prism_rng *rng) {
    return (prism_normal_source){rng, sizeof(*rng), fill_normal};
}
