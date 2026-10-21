/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include <math.h>
#include <prism.h>
#include <stdio.h>
#include <string.h>
int main(int argc, char **argv) {
    enum { N = 257, J = 3 };
    size_t edges[] = {N, 96, 24, 4};
    double times[] = {0.5, 8, 128};
    if (argc > 2 || (argc == 2 && strcmp(argv[1], "fixed") != 0 && strcmp(argv[1], "heat") != 0)) {
        fprintf(stderr, "Usage: %s [fixed|heat]\n", argv[0]);
        return 2;
    }
    prism_mode mode =
        argc == 2 && strcmp(argv[1], "fixed") == 0 ? PRISM_BANDS_FIXED : PRISM_BANDS_HEAT_TELESCOPE;
    prism_config c = {N, J, 1.0, mode, edges, times};
    prism_plan *p = NULL;
    prism_workspace *w = NULL;
    double x[N], components[J * N], residual[N], reconstructed[N];
    for (size_t i = 0; i < N; ++i)
        x[i] = 2 + sin(0.4 * (double)i) + 0.3 * cos(1.7 * (double)i);
    prism_result r = {components, J * N, residual, N, 0, 0, 0, 0};
    int status = prism_plan_create(&c, NULL, NULL, &p);
    if (status == PRISM_OK)
        status = prism_workspace_create(p, &w);
    if (status == PRISM_OK)
        status = prism_decompose(p, w, x, &r);
    if (status == PRISM_OK)
        status = prism_reconstruct(&r, reconstructed, N);
    if (status == PRISM_OK)
        printf("PRISM %s: %zu components + residual, relative error %.3e, max error %.3e\n",
               mode == PRISM_BANDS_FIXED ? "fixed" : "heat", r.count,
               r.reconstruction_relative_error, r.reconstruction_max_abs_error);
    else
        fprintf(stderr, "PRISM: %s\n", prism_status_string(status));
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    return status == PRISM_OK ? 0 : 1;
}
