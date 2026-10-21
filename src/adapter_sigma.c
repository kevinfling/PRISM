/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling
 * SIGMA is a separately supplied optional MIT dependency; see docs/SIGMA.md. */
#include "internal.h"
#include <prism/sigma.h>
#include <string.h>
#ifdef PRISM_HAVE_SIGMA
#include <sigma/rng.h>
#endif

int prism_sigma_available(void) {
#ifdef PRISM_HAVE_SIGMA
    return 1;
#else
    return 0;
#endif
}
const char *prism_sigma_revision(void) {
#ifdef PRISM_HAVE_SIGMA
    return PRISM_SIGMA_REVISION;
#else
    return NULL;
#endif
}
#ifdef PRISM_HAVE_SIGMA
static int valid(const prism_sigma_rng *rng) {
    return rng && (rng->words[0] | rng->words[1] | rng->words[2] | rng->words[3]) &&
           rng->has_normal_spare <= 1u && (!rng->has_normal_spare || isfinite(rng->normal_spare));
}
static void load(const prism_sigma_rng *rng, sigma_state_t *state) {
    memcpy(state->xoshiro.s, rng->words, sizeof(rng->words));
    state->vtable = &sigma_vtable_xoshiro256ss;
    state->spare_normal = rng->normal_spare;
    state->has_spare = rng->has_normal_spare != 0;
}
static void save(prism_sigma_rng *rng, const sigma_state_t *state) {
    memcpy(rng->words, state->xoshiro.s, sizeof(rng->words));
    rng->normal_spare = state->spare_normal;
    rng->has_normal_spare = state->has_spare ? 1u : 0u;
}
#endif
int prism_sigma_seed(prism_sigma_rng *rng, uint64_t seed, uint64_t stream) {
    if (!rng)
        return PRISM_EINVAL;
#ifdef PRISM_HAVE_SIGMA
    sigma_splitmix64_t mixer;
    sigma_splitmix64_seed(&mixer, stream);
    sigma_state_t state;
    sigma_xoshiro256ss_seed(&state.xoshiro, seed ^ sigma_splitmix64_next(&mixer));
    state.spare_normal = 0;
    state.has_spare = false;
    save(rng, &state);
    return PRISM_OK;
#else
    (void)seed;
    (void)stream;
    return PRISM_EBACKEND;
#endif
}
int prism_sigma_next_u64(prism_sigma_rng *rng, uint64_t *out) {
    if (!rng || !out || prism_overlap(rng, sizeof(*rng), out, sizeof(*out)))
        return PRISM_EINVAL;
#ifdef PRISM_HAVE_SIGMA
    if (!valid(rng))
        return PRISM_EINVAL;
    sigma_state_t state;
    load(rng, &state);
    *out = sigma_rng_next_u64(&state);
    save(rng, &state);
    return PRISM_OK;
#else
    return PRISM_EBACKEND;
#endif
}
static int fill_normal(void *opaque, double *out, size_t count, size_t capacity) {
    prism_sigma_rng *rng = opaque;
    if (!rng || !out || !count)
        return PRISM_EINVAL;
    if (capacity < count)
        return PRISM_ECAPACITY;
#ifdef PRISM_HAVE_SIGMA
    if (!valid(rng))
        return PRISM_EINVAL;
    sigma_state_t state;
    load(rng, &state);
    for (size_t i = 0; i < count; ++i)
        out[i] = sigma_rng_next_normal(&state);
    save(rng, &state);
    return PRISM_OK;
#else
    return PRISM_EBACKEND;
#endif
}
prism_normal_source prism_sigma_source(prism_sigma_rng *rng) {
    return (prism_normal_source){rng, sizeof(*rng), fill_normal};
}
