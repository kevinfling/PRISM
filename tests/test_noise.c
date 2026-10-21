/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "check.h"
#include <float.h>
#include <math.h>
#include <prism/noise.h>
#include <stdint.h>
#include <string.h>
#define PI 3.141592653589793238462643383279502884L

static void rng_checks(void) {
    static const uint32_t expected[] = {0xa15c02b7u, 0x7b47f409u, 0xba1d3330u,
                                        0x83d2f293u, 0xbfa4784bu, 0xcbed606eu};
    prism_rng a, b;
    CHECK(prism_rng_seed(&a, 42, 54) == PRISM_OK);
    for (size_t i = 0; i < sizeof(expected) / sizeof(*expected); ++i) {
        uint32_t word;
        CHECK(prism_rng_next_u32(&a, &word) == PRISM_OK);
        CHECK(word == expected[i]);
    }
    CHECK(prism_rng_seed(&a, 713, 51) == PRISM_OK);
    b = a;
    double all[101], chunks[101];
    CHECK(prism_rng_normal(&a, all, 101, 101) == PRISM_OK);
    CHECK(prism_rng_normal(&b, chunks, 3, 101) == PRISM_OK);
    CHECK(prism_rng_normal(&b, chunks + 3, 98, 98) == PRISM_OK);
    CHECK(memcmp(all, chunks, sizeof(all)) == 0);
    CHECK(a.state == b.state && a.has_normal_spare == b.has_normal_spare);
    b = a;
    CHECK(prism_rng_normal(&a, all, 101, 100) == PRISM_ECAPACITY);
    CHECK(a.state == b.state && a.has_normal_spare == b.has_normal_spare);
    CHECK(prism_rng_normal(&a, (double *)(void *)&a, 1, 1) == PRISM_EINVAL);
    CHECK(prism_rng_seed(&a, 0, UINT64_MAX) == PRISM_EINVAL);
    CHECK(prism_rng_seed(&a, 0, 0) == PRISM_OK);
    CHECK(prism_rng_uniform(&a, all, 101, 101) == PRISM_OK);
    for (size_t i = 0; i < 101; ++i)
        CHECK(all[i] >= 0 && all[i] < 1);
    double sum = 0, squares = 0;
    for (size_t j = 0; j < 1000; ++j) {
        CHECK(prism_rng_normal(&a, all, 101, 101) == PRISM_OK);
        for (size_t i = 0; i < 101; ++i) {
            sum += all[i];
            squares += all[i] * all[i];
        }
    }
    CHECK(fabs(sum / 101000) < 0.025);
    CHECK(fabs(squares / 101000 - 1) < 0.04);
}
static void covariance_checks(void) {
    double lags[1024];
    CHECK(prism_fgn_covariance(1024, 0.5, 2, lags, 1024) == PRISM_OK);
    CHECK(lags[0] == 2);
    for (size_t k = 1; k < 1024; ++k)
        CHECK(lags[k] == 0);
    const double hs[] = {0.1, 0.3, 0.7, 0.9};
    for (size_t h = 0; h < 4; ++h) {
        CHECK(prism_fgn_covariance(1024, hs[h], 1.7, lags, 1024) == PRISM_OK);
        for (size_t k = 1; k < 32; ++k) {
            long double p = 2 * (long double)hs[h];
            long double expected =
                0.85L * (powl((long double)k + 1, p) - 2 * powl((long double)k, p) +
                         powl((long double)k - 1, p));
            CHECK(fabsl(lags[k] - expected) < 2e-13L);
        }
        for (size_t k = 1; k < 1024; ++k)
            CHECK(hs[h] < 0.5 ? lags[k] < 0 : lags[k] > 0);
    }
    CHECK(prism_fgn_covariance(1024, nextafter(0.5, 1), 1, lags, 1024) == PRISM_OK);
    CHECK(lags[1023] > 0);
    CHECK(prism_fgn_covariance(1024, nextafter(0.5, 0), 1, lags, 1024) == PRISM_OK);
    CHECK(lags[1023] < 0);
    lags[0] = 27;
    CHECK(prism_fgn_covariance(2, 0.7, 1, lags, 1) == PRISM_ECAPACITY);
    CHECK(lags[0] == 27);
    CHECK(prism_fgn_covariance(SIZE_MAX, 0.7, 1, lags, SIZE_MAX) == PRISM_EOVERFLOW);
    CHECK(prism_fgn_covariance(1, 0, 1, lags, 1) == PRISM_EINVAL);
    CHECK(prism_fgn_covariance(1, 1, 1, lags, 1) == PRISM_EINVAL);
    CHECK(prism_fgn_covariance(1, NAN, 1, lags, 1) == PRISM_EINVAL);
    CHECK(prism_fgn_covariance(1, 0.7, -1, lags, 1) == PRISM_EINVAL);
    CHECK(prism_fgn_covariance(1, 0.7, INFINITY, lags, 1) == PRISM_EINVAL);
}

/* Independent full cosine-matrix square root of the circulant embedding,
 * applied directly in sample space (not a forward/filter/inverse path). */
static void circulant_oracle(size_t n, const double *lags, const double *white, double *out) {
    size_t m = n == 1 ? 1 : 2 * (n - 1);
    long double eigen[64], impulse[64];
    CHECK(m <= 64);
    for (size_t k = 0; k < m; ++k) {
        eigen[k] = 0;
        for (size_t j = 0; j < m; ++j)
            eigen[k] += lags[j < n ? j : m - j] *
                        cosl(2 * PI * (long double)j * (long double)k / (long double)m);
        CHECK(eigen[k] >= 0);
        eigen[k] = sqrtl(eigen[k]);
    }
    for (size_t i = 0; i < m; ++i) {
        impulse[i] = 0;
        for (size_t k = 0; k < m; ++k)
            impulse[i] +=
                eigen[k] * cosl(2 * PI * (long double)i * (long double)k / (long double)m);
        impulse[i] /= (long double)m;
    }
    for (size_t i = 0; i < n; ++i) {
        long double sum = 0;
        for (size_t j = 0; j < m; ++j)
            sum += impulse[(i + m - j) % m] * white[j];
        out[i] = (double)sum;
    }
}
static void oracle_checks(prism_noise_method method) {
    const size_t sizes[] = {1, 2, 3, 8, 17};
    for (size_t s = 0; s < sizeof(sizes) / sizeof(*sizes); ++s) {
        size_t n = sizes[s], m = n == 1 ? 1 : 2 * (n - 1);
        double lags[17], white[32], expected[17], actual[17];
        for (size_t i = 0; i < n; ++i)
            lags[i] = 2.3 * pow(0.6, (double)i);
        prism_noise_config config = {n, method, 0};
        prism_noise_diagnostics diag;
        prism_noise_plan *p = NULL;
        prism_noise_workspace *w = NULL;
        CHECK(prism_noise_plan_create(&config, lags, n, NULL, &p, &diag) == PRISM_OK);
        CHECK(prism_noise_workspace_create(p, &w) == PRISM_OK);
        CHECK(diag.embedding_size == m && diag.rounded_eigenvalue_count == 0);
        for (uint64_t seed = 0; seed < 4; ++seed) {
            prism_rng a, b;
            CHECK(prism_rng_seed(&a, seed, 3) == PRISM_OK);
            b = a;
            CHECK(prism_rng_normal(&b, white, m, m) == PRISM_OK);
            circulant_oracle(n, lags, white, expected);
            CHECK(prism_noise_sample(p, w, &a, actual, n) == PRISM_OK);
            for (size_t i = 0; i < n; ++i)
                CHECK(fabs(actual[i] - expected[i]) < 3e-13);
            CHECK(a.state == b.state);
        }
        prism_noise_workspace_destroy(w);
        prism_noise_plan_destroy(p);
    }
}
static void failure_checks(prism_noise_method method) {
    prism_noise_config config = {3, PRISM_NOISE_CHOLESKY, 0};
    prism_noise_diagnostics diag;
    prism_noise_plan *p = NULL;
    double lags[] = {1.05, 0.9, 0.62};
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_OK);
    prism_noise_plan_destroy(p);
    config.method = method;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_ECOVARIANCE);
    CHECK(!p && fabs(diag.minimum_embedding_eigenvalue + 0.13) < 1e-14);
    lags[0] = 1;
    lags[1] = nextafter(0.5, 1);
    lags[2] = 0;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_ECOVARIANCE);
    config.embedding_relative_tolerance = 1e-12;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_OK);
    CHECK(diag.rounded_eigenvalue_count == 1);
    CHECK(diag.maximum_eigenvalue_adjustment > 0);
    double copied[3];
    CHECK(prism_noise_target_covariance(p, copied, 3) == PRISM_OK);
    CHECK(memcmp(copied, lags, sizeof(lags)) == 0);
    prism_noise_plan_destroy(p);
    config.embedding_relative_tolerance = 0;
    CHECK(prism_noise_plan_create(&config, lags, 2, NULL, &p, &diag) == PRISM_ECAPACITY);
    CHECK(!p && diag.n == 3);
    config.n = SIZE_MAX;
    CHECK(prism_noise_plan_create(&config, lags, SIZE_MAX, NULL, &p, &diag) == PRISM_EOVERFLOW);
    config.n = 3;
    lags[0] = NAN;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_EINVAL);
    lags[0] = -1;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_ECOVARIANCE);
    lags[0] = 1;
    lags[1] = 0.9;
    lags[2] = -0.9;
    config.method = PRISM_NOISE_CHOLESKY;
    CHECK(prism_noise_plan_create(&config, lags, 3, NULL, &p, &diag) == PRISM_ECOVARIANCE);
}

typedef struct {
    size_t calls, live, fail_at;
    int locked;
} allocation_state;
static void *tracked_allocate(void *user, size_t bytes) {
    allocation_state *a = user;
    CHECK(!a->locked);
    if (++a->calls == a->fail_at)
        return NULL;
    void *p = malloc(bytes);
    if (p)
        ++a->live;
    return p;
}
static void tracked_free(void *user, void *p) {
    allocation_state *a = user;
    if (p) {
        CHECK(a->live > 0);
        --a->live;
        free(p);
    }
}
static void allocation_checks(prism_noise_method method) {
    double lags[] = {1, 0.3, 0.09, 0.027};
    prism_noise_config config = {4, method, 0};
    for (size_t fail = 1; fail < 16; ++fail) {
        allocation_state state = {0, 0, fail, 0};
        prism_allocator allocator = {&state, tracked_allocate, tracked_free};
        prism_noise_plan *p = NULL;
        prism_noise_workspace *w = NULL;
        prism_noise_diagnostics diag;
        int status = prism_noise_plan_create(&config, lags, 4, &allocator, &p, &diag);
        if (status == PRISM_OK) {
            status = prism_noise_workspace_create(p, &w);
            if (status == PRISM_OK) {
                double out[4] = {23, 23, 23, 23};
                prism_rng rng, before;
                CHECK(prism_rng_seed(&rng, 42, 9) == PRISM_OK);
                before = rng;
                state.locked = 1;
                CHECK(prism_noise_sample(p, w, &rng, out, 3) == PRISM_ECAPACITY);
                CHECK(out[0] == 23 && rng.state == before.state);
                CHECK(prism_noise_sample(p, w, &rng, out, 4) == PRISM_OK);
                state.locked = 0;
                prism_noise_workspace_destroy(w);
            } else {
                CHECK(status == PRISM_ENOMEM && !w);
            }
            prism_noise_plan_destroy(p);
        } else {
            CHECK(status == PRISM_ENOMEM && !p);
        }
        CHECK(state.live == 0);
    }
}
static void degenerate_checks(prism_noise_method method) {
    double lags[] = {0, 0, 0, 0}, out[4];
    prism_noise_config config = {4, method, 0};
    for (size_t kind = 0; kind < (method == PRISM_NOISE_CHOLESKY ? 2u : 1u); ++kind) {
        if (kind)
            for (size_t i = 0; i < 4; ++i)
                lags[i] = 1;
        prism_noise_plan *p;
        prism_noise_workspace *w;
        prism_noise_diagnostics diag;
        prism_rng rng;
        CHECK(prism_noise_plan_create(&config, lags, 4, NULL, &p, &diag) == PRISM_OK);
        CHECK(prism_noise_workspace_create(p, &w) == PRISM_OK);
        CHECK(prism_rng_seed(&rng, 123, 4) == PRISM_OK);
        CHECK(prism_noise_sample(p, w, &rng, out, 4) == PRISM_OK);
        for (size_t i = 1; i < 4; ++i)
            CHECK(out[i] == out[0]);
        if (!kind)
            CHECK(out[0] == 0);
        prism_noise_workspace_destroy(w);
        prism_noise_plan_destroy(p);
    }
}
static void ownership_checks(prism_noise_method method) {
    double lags[] = {1, 0.3, 0.09, 0.027}, copied[4], out[6];
    prism_noise_config config = {4, method, 0};
    prism_noise_diagnostics diag;
    prism_noise_plan *p, *other;
    prism_noise_workspace *w;
    CHECK(prism_noise_plan_create(&config, lags, 4, NULL, &p, &diag) == PRISM_OK);
    CHECK(prism_noise_plan_create(&config, lags, 4, NULL, &other, &diag) == PRISM_OK);
    CHECK(prism_noise_workspace_create(p, &w) == PRISM_OK);
    lags[0] = 7;
    CHECK(prism_noise_target_covariance(p, copied, 4) == PRISM_OK);
    CHECK(copied[0] == 1);
    prism_rng rng, before;
    CHECK(prism_rng_seed(&rng, 42, 9) == PRISM_OK);
    before = rng;
    for (size_t i = 0; i < 6; ++i)
        out[i] = 23;
    CHECK(prism_noise_sample(other, w, &rng, out, 6) == PRISM_EINVAL);
    CHECK(out[0] == 23 && rng.state == before.state);
    CHECK(prism_noise_sample(p, w, &rng, (double *)(void *)&rng, 4) == PRISM_EINVAL);
    CHECK(rng.state == before.state);
    CHECK(prism_noise_target_covariance(p, out, 3) == PRISM_ECAPACITY);
    CHECK(out[0] == 23);
    prism_rng invalid = {0};
    CHECK(prism_noise_sample(p, w, &invalid, out, 6) == PRISM_EINVAL);
    CHECK(out[0] == 23);
    CHECK(prism_noise_sample(p, w, &rng, out, 6) == PRISM_OK);
    CHECK(out[4] == 23 && out[5] == 23);
    prism_noise_workspace_destroy(w);
    prism_noise_plan_destroy(other);
    prism_noise_plan_destroy(p);
}

/* Known-zero-mean second moment; each selected entry has Gaussian variance
 * (K_ii*K_jj + K_ij^2)/draws. Six-standard-error bounds, fixed seeds. */
static void monte_carlo(prism_noise_method method, size_t n) {
    const size_t draws = 8192;
    const double hs[] = {0.1, 0.3, 0.5, 0.7, 0.9};
    double lags[128], out[128];
    CHECK(n <= 128);
    for (size_t h = 0; h < 5; ++h) {
        CHECK(prism_fgn_covariance(n, hs[h], 1, lags, n) == PRISM_OK);
        prism_noise_config config = {n, method, 0};
        prism_noise_plan *p;
        prism_noise_workspace *w;
        prism_noise_diagnostics diag;
        CHECK(prism_noise_plan_create(&config, lags, n, NULL, &p, &diag) == PRISM_OK);
        CHECK(prism_noise_workspace_create(p, &w) == PRISM_OK);
        prism_rng rng;
        CHECK(prism_rng_seed(&rng, UINT64_C(8675309) + h, (uint64_t)n) == PRISM_OK);
        const size_t left[] = {0, 0, 0, 1, 3, 7};
        size_t right[] = {0, 1, n - 1, 2, n / 2, n - 2};
        double sums[6] = {0}, mean = 0;
        double second_moment[128][128] = {{0}};
        for (size_t s = 0; s < draws; ++s) {
            CHECK(prism_noise_sample(p, w, &rng, out, n) == PRISM_OK);
            mean += out[0];
            for (size_t i = 0; i < n; ++i)
                for (size_t j = 0; j <= i; ++j)
                    second_moment[i][j] += out[i] * out[j];
            for (size_t j = 0; j < 6; ++j)
                sums[j] += out[left[j]] * out[right[j]];
        }
        CHECK(fabs(mean / (double)draws) < 6 / sqrt((double)draws));
        double max_z = 0;
        for (size_t j = 0; j < 6; ++j) {
            size_t lag = right[j] - left[j];
            double se = sqrt((1 + lags[lag] * lags[lag]) / (double)draws);
            double z = fabs(sums[j] / (double)draws - lags[lag]) / se;
            max_z = fmax(max_z, z);
            CHECK(z < 6);
        }
        double squared_defect = 0, squared_target = 0;
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j <= i; ++j) {
                double difference = second_moment[i][j] / (double)draws - lags[i - j];
                double multiplicity = i == j ? 1 : 2;
                squared_defect += multiplicity * difference * difference;
                squared_target += multiplicity * lags[i - j] * lags[i - j];
            }
        double expected_squared_defect = (squared_target + (double)n * (double)n) / (double)draws;
        CHECK(squared_defect < 8 * expected_squared_defect);
        printf("covariance method=%d n=%zu H=%.1f draws=%zu max_selected_entry_z=%.4f "
               "relative_Frobenius=%.6f predicted_rms_relative_Frobenius=%.6f\n",
               (int)method, n, hs[h], draws, max_z, sqrt(squared_defect / squared_target),
               sqrt(expected_squared_defect / squared_target));
        prism_noise_workspace_destroy(w);
        prism_noise_plan_destroy(p);
    }
}
int main(void) {
    rng_checks();
    covariance_checks();
    failure_checks(PRISM_NOISE_CIRCULANT_REFERENCE);
    oracle_checks(PRISM_NOISE_CIRCULANT_REFERENCE);
    allocation_checks(PRISM_NOISE_CHOLESKY);
    allocation_checks(PRISM_NOISE_CIRCULANT_REFERENCE);
    ownership_checks(PRISM_NOISE_CHOLESKY);
    ownership_checks(PRISM_NOISE_CIRCULANT_REFERENCE);
    degenerate_checks(PRISM_NOISE_CHOLESKY);
    degenerate_checks(PRISM_NOISE_CIRCULANT_REFERENCE);
    monte_carlo(PRISM_NOISE_CHOLESKY, 32);
    if (prism_noise_method_available(PRISM_NOISE_CIRCULANT_FFTW)) {
        failure_checks(PRISM_NOISE_CIRCULANT_FFTW);
        oracle_checks(PRISM_NOISE_CIRCULANT_FFTW);
        allocation_checks(PRISM_NOISE_CIRCULANT_FFTW);
        ownership_checks(PRISM_NOISE_CIRCULANT_FFTW);
        degenerate_checks(PRISM_NOISE_CIRCULANT_FFTW);
        monte_carlo(PRISM_NOISE_CIRCULANT_FFTW, 128);
    } else {
        prism_noise_config config = {1, PRISM_NOISE_CIRCULANT_FFTW, 0};
        prism_noise_diagnostics diag;
        prism_noise_plan *p;
        double lag = 1;
        CHECK(prism_noise_plan_create(&config, &lag, 1, NULL, &p, &diag) == PRISM_EBACKEND);
    }
    return 0;
}
