/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "check.h"
#include <float.h>
#include <math.h>
#include <prism/noise.h>
#include <prism/sigma.h>
#include <stdint.h>
#include <string.h>

typedef struct {
    size_t calls;
    int status;
    double value;
} fixture;
static int fill(void *opaque, double *out, size_t count, size_t capacity) {
    fixture *state = opaque;
    ++state->calls;
    CHECK(capacity >= count);
    for (size_t i = 0; i < count; ++i)
        out[i] = state->value;
    return state->status;
}
static void checked_source(void) {
    fixture state = {0, PRISM_OK, 1};
    prism_normal_source source = {&state, sizeof(state), fill};
    double out[8] = {23};
    CHECK(prism_normal_fill(&source, out, 8, 7) == PRISM_ECAPACITY);
    CHECK(state.calls == 0 && out[0] == 23);
    CHECK(prism_normal_fill(&source, out, SIZE_MAX, SIZE_MAX) == PRISM_EOVERFLOW);
    CHECK(state.calls == 0);
    CHECK(prism_normal_fill(&source, (double *)(void *)&state, 1, 1) == PRISM_EINVAL);
    CHECK(prism_normal_fill(&source, (double *)(void *)&source, 1, 1) == PRISM_EINVAL);
    CHECK(state.calls == 0);
    CHECK(prism_normal_fill(&source, out, 8, 8) == PRISM_OK);
    CHECK(state.calls == 1 && out[7] == 1);
    state.value = NAN;
    CHECK(prism_normal_fill(&source, out, 8, 8) == PRISM_ENUMERIC);
    state.status = PRISM_EBACKEND;
    CHECK(prism_normal_fill(&source, out, 8, 8) == PRISM_EBACKEND);
    source.fill_normal = NULL;
    CHECK(prism_normal_fill(&source, out, 8, 8) == PRISM_EINVAL);
}
static void covariance_source(prism_noise_method method) {
    double lags[] = {1, 0.3, 0.09, 0.027};
    prism_noise_config config = {4, method, 0};
    prism_noise_diagnostics diag;
    prism_noise_plan *p;
    prism_noise_workspace *w;
    CHECK(prism_noise_plan_create(&config, lags, 4, NULL, &p, &diag) == PRISM_OK);
    CHECK(prism_noise_workspace_create(p, &w) == PRISM_OK);
    size_t count = prism_noise_innovation_count(p);
    CHECK(count == (method == PRISM_NOISE_CHOLESKY ? 4u : 6u));
    CHECK(prism_noise_innovation_count(NULL) == 0);
    prism_rng rng, clone;
    CHECK(prism_rng_seed(&rng, 42, 54) == PRISM_OK);
    clone = rng;
    double white[6], actual[4], expected[4], old[4];
    prism_normal_source source = prism_rng_source(&rng);
    CHECK(prism_rng_normal(&clone, white, count, 6) == PRISM_OK);
    CHECK(prism_noise_transform(p, w, white, count, expected, 4) == PRISM_OK);
    CHECK(prism_noise_sample_from(p, w, &source, actual, 4) == PRISM_OK);
    CHECK(memcmp(actual, expected, sizeof(actual)) == 0);
    CHECK(prism_rng_seed(&clone, 42, 54) == PRISM_OK);
    CHECK(prism_noise_sample(p, w, &clone, old, 4) == PRISM_OK);
    CHECK(memcmp(actual, old, sizeof(old)) == 0);
    CHECK(prism_noise_transform(p, w, white, count - 1, actual, 4) == PRISM_ECAPACITY);
    CHECK(prism_noise_transform(p, w, white, count, white, 4) == PRISM_EINVAL);
    white[count - 1] = NAN;
    CHECK(prism_noise_transform(p, w, white, count, actual, 4) == PRISM_ENUMERIC);
    CHECK(memcmp(actual, old, sizeof(old)) == 0);
    for (size_t i = 0; i < count; ++i)
        white[i] = DBL_MAX;
    CHECK(prism_noise_transform(p, w, white, count, actual, 4) == PRISM_ENUMERIC);
    CHECK(memcmp(actual, old, sizeof(old)) == 0);
    fixture state = {0, PRISM_ENOMEM, 5};
    source = (prism_normal_source){&state, sizeof(state), fill};
    CHECK(prism_noise_sample_from(p, w, &source, actual, 3) == PRISM_ECAPACITY);
    CHECK(state.calls == 0);
    CHECK(prism_noise_sample_from(p, w, &source, actual, 4) == PRISM_ENOMEM);
    CHECK(state.calls == 1 && memcmp(actual, old, sizeof(old)) == 0);
    state.status = PRISM_OK;
    state.value = INFINITY;
    CHECK(prism_noise_sample_from(p, w, &source, actual, 4) == PRISM_ENUMERIC);
    CHECK(memcmp(actual, old, sizeof(old)) == 0);
    prism_noise_workspace_destroy(w);
    prism_noise_plan_destroy(p);
}
/* Separate specification of SplitMix64 seeding and xoshiro256** words. */
static uint64_t split(uint64_t *x) {
    *x += UINT64_C(0x9e3779b97f4a7c15);
    uint64_t y = *x;
    y = (y ^ (y >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    y = (y ^ (y >> 27)) * UINT64_C(0x94d049bb133111eb);
    return y ^ (y >> 31);
}
static uint64_t rotate(uint64_t x, unsigned n) { return (x << n) | (x >> (64u - n)); }
static uint64_t word(uint64_t s[4]) {
    uint64_t result = rotate(s[1] * 5u, 7) * 9u;
    uint64_t next[] = {s[0] ^ s[1] ^ s[3], s[1] ^ s[0] ^ s[2], s[2] ^ s[0] ^ (s[1] << 17),
                       rotate(s[3] ^ s[1], 45)};
    memcpy(s, next, sizeof(next));
    return result;
}
static void sigma_checks(void) {
    prism_sigma_rng rng;
    double output[101];
    if (!prism_sigma_available()) {
        memset(&rng, 0, sizeof(rng));
        CHECK(prism_sigma_revision() == NULL);
        CHECK(prism_sigma_seed(&rng, 42, 54) == PRISM_EBACKEND);
        uint64_t untouched = 23;
        CHECK(prism_sigma_next_u64(&rng, &untouched) == PRISM_EBACKEND);
        CHECK(untouched == 23);
        prism_normal_source source = prism_sigma_source(&rng);
        output[0] = 23;
        CHECK(prism_normal_fill(&source, output, 1, 101) == PRISM_EBACKEND);
        CHECK(output[0] == 23);
        return;
    }
    CHECK(strcmp(prism_sigma_revision(), "6e9a1ef526365474738d80d4eb71fd945a13abd6") == 0);
    uint64_t mixer = 54, effective = UINT64_C(42) ^ split(&mixer), state[4];
    for (size_t i = 0; i < 4; ++i)
        state[i] = split(&effective);
    CHECK(prism_sigma_seed(&rng, 42, 54) == PRISM_OK);
    for (size_t i = 0; i < 100; ++i) {
        uint64_t actual;
        CHECK(prism_sigma_next_u64(&rng, &actual) == PRISM_OK);
        CHECK(actual == word(state));
    }
    CHECK(prism_sigma_seed(&rng, 42, 54) == PRISM_OK);
    prism_sigma_rng clone = rng;
    prism_normal_source a = prism_sigma_source(&rng), b = prism_sigma_source(&clone);
    double chunks[101];
    CHECK(prism_normal_fill(&a, output, 101, 101) == PRISM_OK);
    CHECK(prism_normal_fill(&b, chunks, 1, 101) == PRISM_OK);
    CHECK(prism_normal_fill(&b, chunks + 1, 2, 100) == PRISM_OK);
    CHECK(prism_normal_fill(&b, chunks + 3, 98, 98) == PRISM_OK);
    CHECK(memcmp(output, chunks, sizeof(output)) == 0);
    CHECK(memcmp(rng.words, clone.words, sizeof(rng.words)) == 0);
    CHECK(rng.has_normal_spare == clone.has_normal_spare);
    clone = rng;
    CHECK(prism_normal_fill(&a, output, 101, 100) == PRISM_ECAPACITY);
    CHECK(memcmp(rng.words, clone.words, sizeof(rng.words)) == 0);
    CHECK(prism_normal_fill(&a, (double *)(void *)&rng, 1, 1) == PRISM_EINVAL);
    /* Empirical moments and a fixed-grid Gaussian CDF check: statistical smoke,
     * not a comprehensive PRNG test or proof of independent normal variates. */
    long double sum = 0, sum2 = 0, sum3 = 0, sum4 = 0;
    size_t below[9] = {0};
    const size_t draws = 202000;
    for (size_t block = 0; block < draws / 101; ++block) {
        CHECK(prism_normal_fill(&a, output, 101, 101) == PRISM_OK);
        for (size_t i = 0; i < 101; ++i) {
            long double v = output[i], v2 = v * v;
            sum += v;
            sum2 += v2;
            sum3 += v2 * v;
            sum4 += v2 * v2;
            for (size_t j = 0; j < 9; ++j)
                if (output[i] <= 0.5 * (double)j - 2)
                    ++below[j];
        }
    }
    CHECK(fabsl(sum / draws) < 6 / sqrt((double)draws));
    CHECK(fabsl(sum2 / draws - 1) < 6 * sqrt(2.0 / (double)draws));
    CHECK(fabsl(sum3 / draws) < 6 * sqrt(15.0 / (double)draws));
    CHECK(fabsl(sum4 / draws - 3) < 6 * sqrt(96.0 / (double)draws));
    double max_z = 0;
    for (size_t j = 0; j < 9; ++j) {
        double threshold = 0.5 * (double)j - 2;
        double p = 0.5 * erfc(-threshold / sqrt(2));
        double z = fabs((double)below[j] / (double)draws - p) / sqrt(p * (1 - p) / (double)draws);
        max_z = fmax(max_z, z);
        CHECK(z < 6);
    }
    printf("SIGMA draws=%zu mean=%.8Lg second=%.8Lg third=%.8Lg fourth=%.8Lg max_cdf_z=%.6g\n",
           draws, sum / draws, sum2 / draws, sum3 / draws, sum4 / draws, max_z);
    /* Exercise the actual optional provider through covariance staging. */
    double lags[4], sample[4];
    prism_noise_config config = {4, PRISM_NOISE_CHOLESKY, 0};
    prism_noise_diagnostics diagnostics;
    prism_noise_plan *plan;
    prism_noise_workspace *workspace;
    CHECK(prism_fgn_covariance(4, 0.7, 1, lags, 4) == PRISM_OK);
    CHECK(prism_noise_plan_create(&config, lags, 4, NULL, &plan, &diagnostics) == PRISM_OK);
    CHECK(prism_noise_workspace_create(plan, &workspace) == PRISM_OK);
    CHECK(prism_noise_sample_from(plan, workspace, &a, sample, 4) == PRISM_OK);
    prism_noise_workspace_destroy(workspace);
    prism_noise_plan_destroy(plan);
}
int main(void) {
    checked_source();
    covariance_source(PRISM_NOISE_CHOLESKY);
    covariance_source(PRISM_NOISE_CIRCULANT_REFERENCE);
    if (prism_noise_method_available(PRISM_NOISE_CIRCULANT_FFTW))
        covariance_source(PRISM_NOISE_CIRCULANT_FFTW);
    sigma_checks();
    return 0;
}
