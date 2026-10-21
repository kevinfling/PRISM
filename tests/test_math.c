/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#include "check.h"
#include <math.h>
#include <prism.h>
#include <string.h>

static void heat_weights(double spacing, const double *times, size_t count, double *weights) {
    prism_config config = {33, count, spacing, PRISM_BANDS_HEAT_TELESCOPE, NULL, times};
    prism_plan *p = NULL;
    CHECK(prism_plan_create(&config, prism_backend_reference(), NULL, &p) == PRISM_OK);
    CHECK(prism_plan_weights(p, weights, (count + 1) * 33) == PRISM_OK);
    prism_plan_destroy(p);
}

static void heat_identities(void) {
    double times[] = {0.3, 2, 9}, scaled_times[] = {4.8, 32, 144};
    double weights[132], scaled[132];
    heat_weights(0.2, times, 3, weights);
    heat_weights(0.8, scaled_times, 3, scaled);
    for (size_t k = 0; k < 33; ++k) {
        long double sum = 0, squares = 0;
        for (size_t j = 0; j < 4; ++j) {
            double w = weights[33 * j + k];
            CHECK(fabs(w - scaled[33 * j + k]) < 1e-14);
            CHECK(w >= 0 && w <= 1);
            sum += w;
            squares += (long double)w * w;
        }
        CHECK(fabsl(sum - 1) < 1e-14L);
        CHECK(squares >= 0.25L - 1e-14L && squares <= 1 + 1e-14L);
    }
    double a[] = {0.7}, b[] = {1.9}, ab[] = {a[0] + b[0]};
    double wa[66], wb[66], wab[66];
    heat_weights(0.2, a, 1, wa);
    heat_weights(0.2, b, 1, wb);
    heat_weights(0.2, ab, 1, wab);
    for (size_t k = 0; k < 33; ++k)
        CHECK(fabs(wa[33 + k] * wb[33 + k] - wab[33 + k]) < 1e-14);
    double close_times[] = {1.0, nextafter(1.0, 2.0)};
    double close_weights[99];
    heat_weights(1, close_times, 2, close_weights);
    /* A tiny positive band must not be lost to subtraction of rounded exponentials. */
    CHECK(close_weights[34] > 0);
    for (size_t k = 0; k < 33; ++k)
        CHECK(fabs(close_weights[k] + close_weights[33 + k] + close_weights[66 + k] - 1) < 1e-14);
}

static void operator_checks(const prism_backend *backend, prism_mode mode) {
    size_t edges[] = {33, 16, 4};
    double times[] = {0.5, 8};
    prism_config config = {33, 2, 1, mode, edges, times};
    prism_plan *p = NULL;
    prism_workspace *w = NULL;
    CHECK(prism_plan_create(&config, backend, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    double x[33], comp[66], residual[33], second[66], second_residual[33];
    prism_result r = {comp, 66, residual, 33, 0, 0, 0, 0};
    prism_result again = {second, 66, second_residual, 33, 0, 0, 0, 0};
    long double input_energy = 0, stacked_energy = 0;
    for (size_t i = 0; i < 33; ++i) {
        x[i] = sin(0.17 * (double)i) + 0.3 * cos(1.2 * (double)i);
        input_energy += (long double)x[i] * x[i];
    }
    CHECK(prism_decompose(p, w, x, &r) == PRISM_OK);
    for (size_t j = 0; j < 3; ++j) {
        const double *component = j == 2 ? residual : comp + 33 * j;
        long double energy = 0;
        for (size_t i = 0; i < 33; ++i)
            energy += (long double)component[i] * component[i];
        CHECK(energy <= input_energy + 1e-12L);
        stacked_energy += energy;
        CHECK(prism_decompose(p, w, component, &again) == PRISM_OK);
        if (mode == PRISM_BANDS_FIXED) {
            const double *projected = j == 2 ? second_residual : second + 33 * j;
            for (size_t i = 0; i < 33; ++i)
                CHECK(fabs(projected[i] - component[i]) < 1e-12);
        }
    }
    CHECK(stacked_energy <= input_energy + 1e-12L);
    CHECK(stacked_energy >= input_energy / 3 - 1e-12L);
    if (mode == PRISM_BANDS_HEAT_TELESCOPE)
        CHECK(stacked_energy < 0.99L * input_energy);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
}

int main(void) {
    heat_identities();
    const prism_backend *backends[] = {prism_backend_reference(), prism_backend_fftw()};
    for (size_t b = 0; b < 2; ++b) {
        if (!backends[b])
            continue;
        operator_checks(backends[b], PRISM_BANDS_FIXED);
        operator_checks(backends[b], PRISM_BANDS_HEAT_TELESCOPE);
    }
    puts("heat semigroup, units, close-time bands, projections, and norm bounds passed");
    return 0;
}
