/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling
 * Fixed protocol; see docs/STOCHASTIC.md. This is a research driver, not a
 * production ensemble decomposition API. CSV on stdout, protocol on stderr. */
#define _POSIX_C_SOURCE 200809L
#include <inttypes.h>
#include <math.h>
#include <prism/noise.h>
#include <prism/sigma.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#define N 64
#define SEEDS 20
#define PI 3.141592653589793238462643383279502884L

static void require(int status) {
    if (status != PRISM_OK) {
        fprintf(stderr, "%s\n", prism_status_string(status));
        exit(1);
    }
}
static void *allocate(void *user, size_t bytes) {
    (void)user;
    return malloc(bytes);
}
static void deallocate(void *user, void *p) {
    (void)user;
    free(p);
}
static double seconds(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t))
        exit(1);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}
static void filter(const prism_backend *backend, void *state, const double *weights,
                   const double *x, double *spectrum, double *out) {
    require(backend->forward(state, x, spectrum));
    for (size_t k = 0; k < N; ++k)
        spectrum[k] *= weights[k];
    require(backend->inverse(state, spectrum, out));
}
static void analytic_moments(const double *weights, const double *lags, double *trace,
                             double *trace_square) {
    /* Build T from a direct orthonormal cosine matrix, then B=T*K*T^T.
     * This calculation does not use either implementation transform. */
    long double q[N][N], t[N][N], tk[N][N], b[N][N];
    for (size_t k = 0; k < N; ++k)
        for (size_t i = 0; i < N; ++i)
            q[k][i] = sqrtl((k ? 2.0L : 1.0L) / N) *
                      cosl(PI * (long double)k * ((long double)i + 0.5L) / N);
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j) {
            t[i][j] = 0;
            for (size_t k = 0; k < N; ++k)
                t[i][j] += q[k][i] * weights[k] * q[k][j];
        }
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j) {
            tk[i][j] = 0;
            for (size_t k = 0; k < N; ++k)
                tk[i][j] += t[i][k] * lags[k > j ? k - j : j - k];
        }
    long double tr = 0, tr2 = 0;
    for (size_t i = 0; i < N; ++i)
        for (size_t j = 0; j < N; ++j) {
            b[i][j] = 0;
            for (size_t k = 0; k < N; ++k)
                b[i][j] += tk[i][k] * t[j][k];
            if (i == j)
                tr += b[i][j];
            tr2 += b[i][j] * b[i][j];
        }
    *trace = (double)tr;
    *trace_square = (double)tr2;
}
int main(int argc, char **argv) {
    if (argc > 3 || (argc >= 2 && strcmp(argv[1], "fftw") && strcmp(argv[1], "reference")) ||
        (argc == 3 && strcmp(argv[2], "pcg") && strcmp(argv[2], "sigma"))) {
        fprintf(stderr, "usage: %s [fftw|reference [pcg|sigma]]\n", argv[0]);
        return 2;
    }
    const prism_backend *backend = prism_backend_fftw();
    prism_noise_method method = PRISM_NOISE_CIRCULANT_FFTW;
    if ((argc >= 2 && !strcmp(argv[1], "reference")) || (argc == 1 && !backend)) {
        backend = prism_backend_reference();
        method = PRISM_NOISE_CHOLESKY;
    }
    if (!backend) {
        fprintf(stderr, "FFTW backend unavailable\n");
        return 2;
    }
    int use_sigma = argc == 3 && !strcmp(argv[2], "sigma");
    if (use_sigma && !prism_sigma_available()) {
        fprintf(stderr, "SIGMA adapter unavailable\n");
        return 2;
    }
    fprintf(stderr, "normal_provider=%s sigma_revision=%s policy=scalar_cached\n",
            use_sigma ? "sigma" : "pcg",
            prism_sigma_revision() ? prism_sigma_revision() : "disabled");
    const double sigma = 0.25, diffusion_time = 16;
    const unsigned q = 3;
    const size_t counts[] = {1, 4, 16, 64, 256};
    const char *controls[] = {"zero", "white", "fgn_h03", "fgn_h07"};
    const double hursts[] = {0.5, 0.5, 0.3, 0.7};
    const uint64_t master_seed = UINT64_C(20261021);
    fprintf(stderr,
            "PRISM %s n=%d dx=1 t=%.17g q=%u sigma=%.17g seeds=%d master_seed=%" PRIu64
            " backend=%s noise_method=%d FFTW=%s planner=ESTIMATE|UNALIGNED\n",
            PRISM_VERSION, N, diffusion_time, q, sigma, SEEDS, master_seed, backend->name,
            (int)method, prism_fftw_version() ? prism_fftw_version() : "disabled");
    fprintf(stderr, "E uses nested prefixes within a replicate. Distinct logical stream IDs; no "
                    "independence proof.\n");
    puts("control,replicate,seed,stream,n,q,t,sigma,E,squared_error,expected_squared_error,"
         "variance_squared_error,linear_identity_max_abs,sampling_seconds,transform_seconds,"
         "aggregation_seconds,setup_seconds,deterministic_seconds");
    prism_allocator allocator = {NULL, allocate, deallocate};
    void *state = NULL;
    double setup_start = seconds();
    require(backend->create(N, &allocator, &state));
    double backend_setup = seconds() - setup_start;
    double x[N], weights[N], spectrum[N], baseline[N];
    for (size_t i = 0; i < N; ++i)
        x[i] = cos((double)(2 * PI * ((long double)i + 0.5L) / N)) +
               0.3 * cos((double)(23 * PI * ((long double)i + 0.5L) / N)) + 0.01 * (double)i;
    for (size_t k = 0; k < N; ++k) {
        double omega = (double)(PI * (long double)k / N);
        weights[k] = pow(-expm1(-diffusion_time * omega * omega), (double)q);
    }
    double start = seconds();
    filter(backend, state, weights, x, spectrum, baseline);
    double deterministic_seconds = seconds() - start;
    for (size_t c = 0; c < 4; ++c) {
        double lags[N];
        require(prism_fgn_covariance(N, hursts[c], c ? 1 : 0, lags, N));
        prism_noise_config config = {N, method, 0};
        prism_noise_plan *plan;
        prism_noise_workspace *workspace;
        prism_noise_diagnostics diag;
        start = seconds();
        require(prism_noise_plan_create(&config, lags, N, NULL, &plan, &diag));
        require(prism_noise_workspace_create(plan, &workspace));
        double setup_seconds = backend_setup + seconds() - start;
        double trace, trace_square;
        analytic_moments(weights, lags, &trace, &trace_square);
        double totals[5] = {0};
        for (size_t seed_index = 0; seed_index < SEEDS; ++seed_index) {
            uint64_t stream = (uint64_t)(4 * seed_index + c);
            prism_rng rng;
            prism_sigma_rng srng;
            prism_normal_source source;
            if (use_sigma) {
                require(prism_sigma_seed(&srng, master_seed, stream));
                source = prism_sigma_source(&srng);
            } else {
                require(prism_rng_seed(&rng, master_seed, stream));
                source = prism_rng_source(&rng);
            }
            double noise[N], trial[N], filtered[N], sum[N] = {0}, noise_sum[N] = {0};
            double sampling_seconds = 0, transform_seconds = 0, aggregation_seconds = 0;
            size_t checkpoint = 0;
            for (size_t e = 1; e <= counts[4]; ++e) {
                start = seconds();
                require(prism_noise_sample_from(plan, workspace, &source, noise, N));
                sampling_seconds += seconds() - start;
                start = seconds();
                for (size_t i = 0; i < N; ++i)
                    trial[i] = x[i] + sigma * noise[i];
                filter(backend, state, weights, trial, spectrum, filtered);
                transform_seconds += seconds() - start;
                start = seconds();
                for (size_t i = 0; i < N; ++i) {
                    sum[i] += filtered[i];
                    noise_sum[i] += noise[i];
                }
                aggregation_seconds += seconds() - start;
                if (e == counts[checkpoint]) {
                    double error = 0, defect = 0;
                    for (size_t i = 0; i < N; ++i)
                        trial[i] = noise_sum[i] / (double)e;
                    filter(backend, state, weights, trial, spectrum, filtered);
                    for (size_t i = 0; i < N; ++i) {
                        double actual = sum[i] / (double)e;
                        double difference = actual - baseline[i];
                        error += difference * difference;
                        defect = fmax(defect, fabs(actual - baseline[i] - sigma * filtered[i]));
                    }
                    if (defect > 2e-12 || !isfinite(error)) {
                        fprintf(stderr, "fixed-linear identity failed\n");
                        return 1;
                    }
                    double scale = sigma * sigma / (double)e;
                    double expected = scale * trace;
                    double variance = 2 * scale * scale * trace_square;
                    totals[checkpoint] += error;
                    printf(
                        "%s,%zu,%" PRIu64 ",%" PRIu64
                        ",%d,%u,%.17g,%.17g,%zu,%.17g,%.17g,%.17g,%.17g,%.9g,%.9g,%.9g,%.9g,%.9g\n",
                        controls[c], seed_index, master_seed, stream, N, q, diffusion_time, sigma,
                        e, error, expected, variance, defect, sampling_seconds, transform_seconds,
                        aggregation_seconds, setup_seconds, deterministic_seconds);
                    ++checkpoint;
                }
            }
        }
        for (size_t j = 0; j < 5; ++j) {
            double scale = sigma * sigma / (double)counts[j];
            double expected = scale * trace;
            double standard_error = sqrt(2 * scale * scale * trace_square / SEEDS);
            double observed = totals[j] / SEEDS;
            if (fabs(observed - expected) > 6 * standard_error + 1e-24) {
                fprintf(stderr,
                        "ensemble second moment failed: %s E=%zu observed=%.17g expected=%.17g "
                        "se=%.17g\n",
                        controls[c], counts[j], observed, expected, standard_error);
                return 1;
            }
        }
        fprintf(stderr, "%s embedding=%zu min_eigenvalue=%.17g rounded=%zu\n", controls[c],
                diag.embedding_size, diag.minimum_embedding_eigenvalue,
                diag.rounded_eigenvalue_count);
        prism_noise_workspace_destroy(workspace);
        prism_noise_plan_destroy(plan);
    }
    backend->destroy(state);
    return 0;
}
