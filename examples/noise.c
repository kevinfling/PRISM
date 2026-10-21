/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include <prism/noise.h>
#include <stdio.h>
int main(void) {
    double lags[32], sample[32];
    prism_noise_config config = {32, PRISM_NOISE_CHOLESKY, 0};
    prism_noise_diagnostics diagnostics;
    prism_noise_plan *plan = NULL;
    prism_noise_workspace *workspace = NULL;
    prism_rng rng;
    int status = prism_fgn_covariance(32, 0.7, 1, lags, 32);
    if (!status)
        status = prism_noise_plan_create(&config, lags, 32, NULL, &plan, &diagnostics);
    if (!status)
        status = prism_noise_workspace_create(plan, &workspace);
    if (!status)
        status = prism_rng_seed(&rng, 42, 54);
    if (!status)
        status = prism_noise_sample(plan, workspace, &rng, sample, 32);
    if (!status)
        for (size_t i = 0; i < 32; ++i)
            printf("%zu,%.17g\n", i, sample[i]);
    else
        fprintf(stderr, "%s\n", prism_status_string(status));
    prism_noise_workspace_destroy(workspace);
    prism_noise_plan_destroy(plan);
    return status ? 1 : 0;
}
