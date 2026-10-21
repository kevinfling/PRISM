/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_RNG_H
#define PRISM_RNG_H
#include <prism/innovation.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Caller-owned PCG-XSH-RR 64/32 state with a Box-Muller spare. Seed before use.
 * Treat fields as private. A generator must not be shared by concurrent calls. */
typedef struct {
    uint64_t state, increment;
    double normal_spare;
    unsigned has_normal_spare;
} prism_rng;

/* Stream identifiers must fit in 63 bits. Equal seed/stream pairs restart the
 * same integer sequence. Stream IDs identify logical trials, not worker IDs;
 * different stream parameters are not a proof of statistical independence. */
int prism_rng_seed(prism_rng *, uint64_t seed, uint64_t stream);
/* Returned descriptor borrows rng; construction does not consume its state. */
prism_normal_source prism_rng_source(prism_rng *rng);
int prism_rng_next_u32(prism_rng *, uint32_t *out);
/* Uniform values are in [0,1). Normals use Box-Muller and libm, so cross-platform
 * bitwise identity is promised for integer outputs only, not normal samples.
 * count must be positive; capacities count doubles. Buffers must not alias state. */
int prism_rng_uniform(prism_rng *, double *out, size_t count, size_t capacity);
int prism_rng_normal(prism_rng *, double *out, size_t count, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
