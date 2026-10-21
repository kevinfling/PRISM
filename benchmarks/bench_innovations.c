/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <prism/noise.h>
#include <prism/sigma.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/utsname.h>
#include <time.h>
#ifndef PRISM_BENCH_BUILD_TYPE
#define PRISM_BENCH_BUILD_TYPE "unknown"
#endif
static double now(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_MONOTONIC, &t))
        exit(1);
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}
static int parse(const char *text, size_t *out) {
    size_t n = 0;
    if (!*text)
        return 0;
    for (const char *p = text; *p; ++p) {
        if (*p < '0' || *p > '9' || n > (SIZE_MAX - (size_t)(*p - '0')) / 10)
            return 0;
        n = 10 * n + (size_t)(*p - '0');
    }
    *out = n;
    return 1;
}
static int seed_source(int sigma, prism_rng *pcg, prism_sigma_rng *srng,
                       prism_normal_source *source) {
    if (sigma) {
        *source = prism_sigma_source(srng);
        return prism_sigma_seed(srng, 20261021, 54);
    }
    *source = prism_rng_source(pcg);
    return prism_rng_seed(pcg, 20261021, 54);
}
int main(int argc, char **argv) {
    size_t n = 257, repetitions = 100;
    if (argc < 2 || argc > 4 || (strcmp(argv[1], "pcg") && strcmp(argv[1], "sigma")) ||
        (argc > 2 && !parse(argv[2], &n)) || (argc > 3 && !parse(argv[3], &repetitions)) || n < 2 ||
        !repetitions || n > SIZE_MAX / sizeof(double) / 3) {
        fprintf(stderr, "usage: %s pcg|sigma [n>=2 [repetitions>=1]]\n", argv[0]);
        return 2;
    }
    int sigma = !strcmp(argv[1], "sigma");
    if (sigma && !prism_sigma_available()) {
        fprintf(stderr, "SIGMA adapter disabled\n");
        return 2;
    }
    const prism_backend *backend = prism_backend_fftw();
    prism_noise_method method = PRISM_NOISE_CIRCULANT_FFTW;
    if (!backend) {
        backend = prism_backend_reference();
        method = PRISM_NOISE_CHOLESKY;
    }
    double *lags = malloc(n * sizeof(double)), *sample = malloc(n * sizeof(double));
    double *signal = malloc(n * sizeof(double)), *components = malloc(3 * n * sizeof(double));
    double *residual = malloc(n * sizeof(double));
    prism_plan *p = NULL;
    prism_workspace *w = NULL;
    prism_noise_plan *np = NULL;
    prism_noise_workspace *nw = NULL;
    int status = PRISM_ENOMEM;
    if (!lags || !sample || !signal || !components || !residual)
        goto done;
    double setup_start = now();
    double times[] = {0.5, 8, 128};
    prism_config config = {n, 3, 1, PRISM_BANDS_HEAT_TELESCOPE, NULL, times};
    prism_noise_config noise_config = {n, method, 0};
    prism_noise_diagnostics diagnostics;
    status = prism_plan_create(&config, backend, NULL, &p);
    if (status)
        goto done;
    status = prism_workspace_create(p, &w);
    if (status)
        goto done;
    status = prism_fgn_covariance(n, 0.7, 1, lags, n);
    if (status)
        goto done;
    status = prism_noise_plan_create(&noise_config, lags, n, NULL, &np, &diagnostics);
    if (status)
        goto done;
    status = prism_noise_workspace_create(np, &nw);
    if (status)
        goto done;
    double setup_seconds = now() - setup_start;
    prism_rng pcg;
    prism_sigma_rng srng;
    prism_normal_source source;
    status = seed_source(sigma, &pcg, &srng, &source);
    if (status)
        goto done;
    for (size_t i = 0; i < 3; ++i) {
        status = prism_normal_fill(&source, sample, n, n);
        if (status)
            goto done;
    }
    double start = now();
    for (size_t i = 0; i < repetitions; ++i) {
        status = prism_normal_fill(&source, sample, n, n);
        if (status)
            goto done;
    }
    double fill_seconds = now() - start;
    long double checksum = 0;
    for (size_t i = 0; i < n; ++i)
        checksum += sample[i];
    status = seed_source(sigma, &pcg, &srng, &source);
    if (status)
        goto done;
    prism_result result = {components, 3 * n, residual, n, 0, 0, 0, 0};
    /* Noise sampling + perturbation + full decomposition including diagnostics.
     * The fixed clean input is precomputed outside the timed region. */
    for (size_t i = 0; i < n; ++i)
        lags[i] = sin(0.4 * (double)i) + 0.001 * (double)i;
    for (size_t i = 0; i < 3; ++i) {
        status = prism_noise_sample_from(np, nw, &source, sample, n);
        if (status)
            goto done;
        for (size_t j = 0; j < n; ++j)
            signal[j] = lags[j] + 0.25 * sample[j];
        status = prism_decompose(p, w, signal, &result);
        if (status)
            goto done;
    }
    start = now();
    for (size_t i = 0; i < repetitions; ++i) {
        status = prism_noise_sample_from(np, nw, &source, sample, n);
        if (status)
            goto done;
        for (size_t j = 0; j < n; ++j)
            signal[j] = lags[j] + 0.25 * sample[j];
        status = prism_decompose(p, w, signal, &result);
        if (status)
            goto done;
    }
    double pipeline_seconds = now() - start;
    struct utsname machine;
    if (!uname(&machine))
        printf("# system=%s %s %s\n", machine.sysname, machine.release, machine.machine);
    printf("# prism=%s build=%s compiler=%s\n", PRISM_VERSION, PRISM_BENCH_BUILD_TYPE, __VERSION__);
    printf("# sigma_revision=%s normal_policy=scalar_cached threads=1 seed=20261021 stream=54 "
           "warmups=3\n",
           prism_sigma_revision() ? prism_sigma_revision() : "disabled");
    printf("# backend=%s noise_method=%d H=0.7 variance=1 sigma=0.25 embedding=%zu rounded=%zu\n",
           backend->name, (int)method, diagnostics.embedding_size,
           diagnostics.rounded_eigenvalue_count);
    puts("provider,n,repetitions,setup_seconds,normal_fill_seconds,normal_ns_per_value,pipeline_"
         "seconds,pipeline_seconds_per_call,last_fill_checksum,reconstruction_relative_error");
    printf("%s,%zu,%zu,%.9g,%.9g,%.9g,%.9g,%.9g,%.17Lg,%.9g\n", argv[1], n, repetitions,
           setup_seconds, fill_seconds, 1e9 * fill_seconds / (double)n / (double)repetitions,
           pipeline_seconds, pipeline_seconds / (double)repetitions, checksum,
           result.reconstruction_relative_error);
done:
    if (status)
        fprintf(stderr, "%s\n", prism_status_string(status));
    prism_noise_workspace_destroy(nw);
    prism_noise_plan_destroy(np);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    free(residual);
    free(components);
    free(signal);
    free(sample);
    free(lags);
    return status ? 1 : 0;
}
