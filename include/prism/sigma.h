/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_SIGMA_H
#define PRISM_SIGMA_H
#include <prism/innovation.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Caller-owned adapter state; fields are private. No SIGMA public header or
 * ABI leaks into PRISM. Seed before use; keep separate state per concurrent call. */
typedef struct {
    uint64_t words[4];
    double normal_spare;
    unsigned has_normal_spare;
} prism_sigma_rng;
int prism_sigma_available(void);
/* Pinned source revision when enabled, NULL otherwise. */
const char *prism_sigma_revision(void);
/* The 64-bit seed and stream map to an effective seed by XOR with one SplitMix64
 * output initialized from stream. Fixed seed + distinct stream gives distinct
 * effective seeds; pairs with different seeds can collide. This is a reproducible
 * indexing policy, not jump-ahead separation or a proof of independence. */
int prism_sigma_seed(prism_sigma_rng *, uint64_t seed, uint64_t stream);
int prism_sigma_next_u64(prism_sigma_rng *, uint64_t *out);
/* Scalar cached normal sequence, including for buffer fills. Deliberately does
 * not select SIGMA SIMD fills or sigma_rng_fill_normal. Chunking is invariant
 * within the same build/libm. Disabled adapters return PRISM_EBACKEND on draws. */
prism_normal_source prism_sigma_source(prism_sigma_rng *);
#ifdef __cplusplus
}
#endif
#endif
