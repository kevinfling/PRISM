/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <prism.h>
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
    if (clock_gettime(CLOCK_MONOTONIC, &t) != 0) {
        perror("clock_gettime");
        exit(1);
    }
    return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

static int parse_size(const char *text, size_t *out) {
    size_t value = 0;
    if (!*text)
        return 0;
    for (const char *p = text; *p; ++p) {
        if (*p < '0' || *p > '9')
            return 0;
        size_t digit = (size_t)(*p - '0');
        if (value > (SIZE_MAX - digit) / 10)
            return 0;
        value = 10 * value + digit;
    }
    *out = value;
    return 1;
}

int main(int argc, char **argv) {
    const prism_backend *b = prism_backend_fftw();
    if (!b)
        b = prism_backend_reference();
    if (argc > 4 ||
        (argc > 1 && strcmp(argv[1], "reference") != 0 && strcmp(argv[1], "fftw") != 0)) {
        fprintf(stderr, "Usage: %s [fftw|reference [n [repetitions]]]\n", argv[0]);
        return 2;
    }
    if (argc > 1)
        b = strcmp(argv[1], "reference") == 0 ? prism_backend_reference() : prism_backend_fftw();
    if (!b) {
        fprintf(stderr, "PRISM: requested FFTW backend is disabled\n");
        return 2;
    }
    size_t n = b == prism_backend_reference() ? 512 : 4096, repetitions = 20;
    if ((argc > 2 && (!parse_size(argv[2], &n) || n < 2)) ||
        (argc > 3 && (!parse_size(argv[3], &repetitions) || !repetitions)) ||
        n > SIZE_MAX / sizeof(double) / 3) {
        fprintf(stderr, "PRISM: n must be >= 2, repetitions >= 1, and sizes representable\n");
        return 2;
    }
    double *x = malloc(n * sizeof(double)), *comp = malloc(3 * n * sizeof(double)),
           *res = malloc(n * sizeof(double));
    prism_plan *p = NULL;
    prism_workspace *w = NULL;
    int status = PRISM_ENOMEM;
    if (!x || !comp || !res)
        goto cleanup;
    for (size_t i = 0; i < n; ++i)
        x[i] = sin(0.4 * (double)i) + 0.001 * (double)i;
    double times[] = {0.5, 8, 128};
    prism_config c = {n, 3, 1, PRISM_BANDS_HEAT_TELESCOPE, NULL, times};
    prism_result r = {comp, 3 * n, res, n, 0, 0, 0, 0};
    double start = now();
    status = prism_plan_create(&c, b, NULL, &p);
    if (status != PRISM_OK)
        goto cleanup;
    status = prism_workspace_create(p, &w);
    if (status != PRISM_OK)
        goto cleanup;
    double setup = now() - start;
    start = now();
    status = prism_decompose(p, w, x, &r);
    if (status != PRISM_OK)
        goto cleanup;
    double first = now() - start;
    for (size_t i = 0; i < 3; ++i) {
        status = prism_decompose(p, w, x, &r);
        if (status != PRISM_OK)
            goto cleanup;
    }
    start = now();
    for (size_t i = 0; i < repetitions; ++i) {
        status = prism_decompose(p, w, x, &r);
        if (status != PRISM_OK)
            goto cleanup;
    }
    double warm = (now() - start) / (double)repetitions;
    struct utsname machine;
    if (uname(&machine) == 0)
        printf("# system=%s %s %s\n", machine.sysname, machine.release, machine.machine);
    FILE *cpu = fopen("/proc/cpuinfo", "r");
    if (cpu) {
        char line[512];
        while (fgets(line, sizeof(line), cpu)) {
            if (strncmp(line, "model name", 10) == 0 || strncmp(line, "Hardware", 8) == 0) {
                printf("# cpu=%s", line);
                break;
            }
        }
        fclose(cpu);
    }
    printf("# prism=%s, build=%s\n# compiler=%s\n# fftw_linked_version=%s\n", PRISM_VERSION,
           PRISM_BENCH_BUILD_TYPE, __VERSION__, prism_fftw_version());
    printf("# planner=%s, execution_threads=1, rng=none, warmup_calls=3\n",
           b == prism_backend_fftw() ? "ESTIMATE|UNALIGNED" : "not_applicable");
    printf("backend,n,components,repetitions,setup_seconds,first_seconds,warm_seconds,relative_"
           "error,max_abs_error\n");
    printf("%s,%zu,3,%zu,%.9g,%.9g,%.9g,%.9g,%.9g\n", b->name, n, repetitions, setup, first, warm,
           r.reconstruction_relative_error, r.reconstruction_max_abs_error);
cleanup:
    if (status != PRISM_OK)
        fprintf(stderr, "PRISM: %s\n", prism_status_string(status));
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    free(res);
    free(comp);
    free(x);
    return status == PRISM_OK ? 0 : 1;
}
