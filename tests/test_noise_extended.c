/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "check.h"
#include <math.h>
#include <prism/noise.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#define MAX_N 513
#define DRAWS 4096

static void dense_fallback(void) {
    /* Keep the dependency-free test practical. Dense Cholesky at the largest
     * awkward size supplements the smaller five-H Monte Carlo suite. */
    double lags[MAX_N], white[MAX_N] = {0}, actual[MAX_N];
    CHECK(prism_fgn_covariance(MAX_N, 0.7, 1, lags, MAX_N) == PRISM_OK);
    prism_noise_config config = {MAX_N, PRISM_NOISE_CHOLESKY, 0};
    prism_noise_diagnostics diagnostics;
    prism_noise_plan *plan;
    prism_noise_workspace *workspace;
    CHECK(prism_noise_plan_create(&config, lags, MAX_N, NULL, &plan, &diagnostics) == PRISM_OK);
    CHECK(prism_noise_workspace_create(plan, &workspace) == PRISM_OK);
    white[0] = 1;
    CHECK(prism_noise_transform(plan, workspace, white, MAX_N, actual, MAX_N) == PRISM_OK);
    /* For unit lag-zero variance, the first Cholesky column equals the first
     * target covariance column exactly up to stored-factor rounding. */
    for (size_t i = 0; i < MAX_N; ++i)
        CHECK(fabs(actual[i] - lags[i]) < 1e-15);
    CHECK(diagnostics.minimum_cholesky_pivot > 0);
    prism_noise_workspace_destroy(workspace);
    prism_noise_plan_destroy(plan);
    puts("dense n=513 H=0.7 first square-root column agrees with target column");
}
static void check_case(size_t n, double h, size_t h_index, FILE *csv) {
    double lags[MAX_N], output[MAX_N], sums[MAX_N] = {0};
    CHECK(prism_fgn_covariance(n, h, 1, lags, MAX_N) == PRISM_OK);
    /* Dense Toeplitz oracle lags from direct powers, not the implementation's
     * stable series; absolute tolerances account for this oracle's cancellation. */
    for (size_t k = 1; k < n; ++k) {
        long double p = 2 * (long double)h, lag = (long double)k;
        long double target = (powl(lag + 1, p) - 2 * powl(lag, p) + powl(lag - 1, p)) / 2;
        CHECK(fabsl(target - lags[k]) < 2e-11L);
    }
    prism_noise_config config = {n, PRISM_NOISE_CIRCULANT_FFTW, 0};
    prism_noise_diagnostics diagnostics;
    prism_noise_plan *plan;
    prism_noise_workspace *workspace;
    CHECK(prism_noise_plan_create(&config, lags, n, NULL, &plan, &diagnostics) == PRISM_OK);
    CHECK(prism_noise_workspace_create(plan, &workspace) == PRISM_OK);
    CHECK(diagnostics.rounded_eigenvalue_count == 0);
    size_t m = prism_noise_innovation_count(plan);
    double *white = calloc(m, sizeof(double)), *factor = malloc(n * m * sizeof(double));
    CHECK(white && factor);
    /* Probe every innovation basis vector and recover the implemented cropped
     * square root R. Compare R*R^T with dense Toeplitz entries, with no random
     * numbers or statistical tolerance in this normalization/factor check. */
    for (size_t k = 0; k < m; ++k) {
        if (k)
            white[k - 1] = 0;
        white[k] = 1;
        CHECK(prism_noise_transform(plan, workspace, white, m, output, n) == PRISM_OK);
        for (size_t i = 0; i < n; ++i)
            factor[i * m + k] = output[i];
    }
    double max_factor_error = 0;
    for (size_t i = 0; i < n; ++i) {
        const size_t js[] = {0, i, n - 1};
        for (size_t pair = 0; pair < 3; ++pair) {
            size_t j = js[pair];
            long double covariance = 0;
            for (size_t k = 0; k < m; ++k)
                covariance += (long double)factor[i * m + k] * factor[j * m + k];
            double defect = (double)fabsl(covariance - lags[i > j ? i - j : j - i]);
            max_factor_error = fmax(max_factor_error, defect);
            CHECK(defect < 3e-12);
        }
    }
    prism_rng rng;
    uint64_t seed = UINT64_C(2026102200) + (uint64_t)h_index;
    CHECK(prism_rng_seed(&rng, seed, (uint64_t)n) == PRISM_OK);
    double mean = 0;
    for (size_t s = 0; s < DRAWS; ++s) {
        CHECK(prism_noise_sample(plan, workspace, &rng, output, n) == PRISM_OK);
        mean += output[0];
        for (size_t k = 0; k < n; ++k)
            sums[k] += output[0] * output[k];
    }
    CHECK(fabs(mean / DRAWS) < 7.5 / sqrt(DRAWS));
    double max_z = 0;
    size_t outside = 0;
    for (size_t k = 0; k < n; ++k) {
        double estimate = sums[k] / DRAWS;
        double se = sqrt((1 + lags[k] * lags[k]) / DRAWS);
        double z = fabs(estimate - lags[k]) / se;
        max_z = fmax(max_z, z);
        if (z > 1.96)
            ++outside;
        CHECK(z < 7.5);
        if (csv)
            fprintf(csv, "%zu,%.1f,%d,%llu,%zu,%zu,%.17g,%.17g,%.17g,%.17g,%.17g\n", n, h, DRAWS,
                    (unsigned long long)seed, n, k, lags[k], estimate, se, estimate - 1.96 * se,
                    estimate + 1.96 * se);
    }
    printf("circulant n=%zu H=%.1f embedding=%zu max_factor_error=%.6g max_lag_z=%.6g "
           "outside_pointwise95=%zu/%zu\n",
           n, h, m, max_factor_error, max_z, outside, n);
    free(factor);
    free(white);
    prism_noise_workspace_destroy(workspace);
    prism_noise_plan_destroy(plan);
}
int main(int argc, char **argv) {
    CHECK(argc <= 2);
    FILE *csv = NULL;
    if (argc == 2) {
        csv = fopen(argv[1], "w");
        CHECK(csv);
        fputs("n,H,draws,seed,stream,lag,target_covariance,estimate,theoretical_se,approx95_lower,"
              "approx95_upper\n",
              csv);
    }
    dense_fallback();
    if (prism_noise_method_available(PRISM_NOISE_CIRCULANT_FFTW)) {
        const size_t sizes[] = {257, 511, 513};
        const double hursts[] = {0.1, 0.3, 0.5, 0.7, 0.9};
        for (size_t s = 0; s < 3; ++s)
            for (size_t h = 0; h < 5; ++h)
                check_case(sizes[s], hursts[h], h, csv);
    }
    if (csv)
        CHECK(fclose(csv) == 0);
    return 0;
}
