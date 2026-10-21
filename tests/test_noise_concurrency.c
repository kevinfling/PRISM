/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "check.h"
#include <prism/noise.h>
#include <prism/sigma.h>
#include <pthread.h>
#include <string.h>
#define N 17
#define JOBS 16
#define THREADS 7

typedef struct {
    const prism_noise_plan *plan;
    prism_noise_workspace *workspace;
    double (*out)[N];
    size_t first, stride;
    int sigma;
} task;
static void *run(void *arg) {
    task *t = arg;
    for (size_t job = t->first; job < JOBS; job += t->stride) {
        /* Every job is one logical sample, independent of worker assignment. */
        size_t sample_id = JOBS - 1 - job;
        prism_rng pcg;
        prism_sigma_rng sigma;
        prism_normal_source source;
        if (t->sigma) {
            CHECK(prism_sigma_seed(&sigma, 20261021, (uint64_t)sample_id) == PRISM_OK);
            source = prism_sigma_source(&sigma);
        } else {
            CHECK(prism_rng_seed(&pcg, 20261021, (uint64_t)sample_id) == PRISM_OK);
            source = prism_rng_source(&pcg);
        }
        CHECK(prism_noise_sample_from(t->plan, t->workspace, &source, t->out[sample_id], N) ==
              PRISM_OK);
    }
    return NULL;
}
static void check_method(prism_noise_method method, int sigma) {
    double lags[N], serial[JOBS][N], parallel[JOBS][N];
    prism_noise_config config = {N, method, 0};
    prism_noise_diagnostics diag;
    prism_noise_plan *p;
    prism_noise_workspace *workspaces[THREADS];
    CHECK(prism_fgn_covariance(N, 0.7, 1, lags, N) == PRISM_OK);
    CHECK(prism_noise_plan_create(&config, lags, N, NULL, &p, &diag) == PRISM_OK);
    for (size_t i = 0; i < THREADS; ++i)
        CHECK(prism_noise_workspace_create(p, &workspaces[i]) == PRISM_OK);
    task single = {p, workspaces[0], serial, 0, 1, sigma};
    (void)run(&single);
    const size_t counts[] = {1, 2, 4, 7};
    for (size_t c = 0; c < 4; ++c) {
        task tasks[THREADS];
        pthread_t threads[THREADS];
        for (size_t i = 0; i < counts[c]; ++i) {
            tasks[i] = (task){p, workspaces[i], parallel, i, counts[c], sigma};
            CHECK(pthread_create(&threads[i], NULL, run, &tasks[i]) == 0);
        }
        for (size_t i = 0; i < counts[c]; ++i)
            CHECK(pthread_join(threads[i], NULL) == 0);
        CHECK(memcmp(serial, parallel, sizeof(serial)) == 0);
    }
    for (size_t i = 0; i < THREADS; ++i)
        prism_noise_workspace_destroy(workspaces[i]);
    prism_noise_plan_destroy(p);
}
int main(void) {
    for (int sigma = 0; sigma <= prism_sigma_available(); ++sigma) {
        check_method(PRISM_NOISE_CHOLESKY, sigma);
        check_method(PRISM_NOISE_CIRCULANT_REFERENCE, sigma);
        if (prism_noise_method_available(PRISM_NOISE_CIRCULANT_FFTW))
            check_method(PRISM_NOISE_CIRCULANT_FFTW, sigma);
    }
    return 0;
}
