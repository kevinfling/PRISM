#include "check.h"
#include <math.h>
#include <prism.h>
#include <stdint.h>
#include <string.h>
static uint64_t seed = UINT64_C(0x27ca1c847e215c93);
static double uniform(void) {
    seed = seed * UINT64_C(6364136223846793005) + UINT64_C(1442695040888963407);
    return (double)(seed >> 11) * 0x1p-53;
}
static void property_test(const prism_backend *b) {
    for (size_t trial = 0; trial < 100; ++trial) {
        size_t n = 2 + (size_t)(uniform() * 62), edges[3] = {n, n / 2, 1};
        size_t count = n > 3 ? 2 : 1;
        double times[2] = {1e-12 + uniform(), 4 + uniform()};
        prism_config c = {n,
                          count,
                          0.001 + 10 * uniform(),
                          trial % 2 ? PRISM_BANDS_FIXED : PRISM_BANDS_HEAT_TELESCOPE,
                          edges,
                          times};
        prism_plan *p = NULL;
        prism_workspace *w = NULL;
        CHECK(prism_plan_create(&c, b, NULL, &p) == PRISM_OK);
        CHECK(prism_workspace_create(p, &w) == PRISM_OK);
        double x[64], y[64], sum[64], cx[128], cy[128], cs[128], rx[64], ry[64], rs[64];
        prism_result a = {cx, 128, rx, 64, 0, 0, 0, 0}, d = {cy, 128, ry, 64, 0, 0, 0, 0},
                     s = {cs, 128, rs, 64, 0, 0, 0, 0};
        for (size_t i = 0; i < n; ++i) {
            x[i] = uniform() - 0.5;
            y[i] = uniform() - 0.5;
            sum[i] = 2 * x[i] - 0.3 * y[i] + 4;
        }
        CHECK(prism_decompose(p, w, x, &a) == PRISM_OK);
        CHECK(prism_decompose(p, w, y, &d) == PRISM_OK);
        CHECK(prism_decompose(p, w, sum, &s) == PRISM_OK);
        CHECK(s.reconstruction_relative_error < 1e-10);
        for (size_t i = 0; i < count * n; ++i)
            CHECK(fabs(cs[i] - (2 * cx[i] - 0.3 * cy[i])) < 1e-12);
        for (size_t i = 0; i < n; ++i)
            CHECK(fabs(rs[i] - (2 * rx[i] - 0.3 * ry[i] + 4)) < 1e-12);
        /* Energy with all cross terms, rather than assuming orthogonality. */
        long double total = 0, diagonal = 0, cross = 0;
        for (size_t i = 0; i < n; ++i) {
            total += (long double)x[i] * x[i];
            diagonal += (long double)rx[i] * rx[i];
            for (size_t j = 0; j < count; ++j) {
                diagonal += (long double)cx[j * n + i] * cx[j * n + i];
                cross += 2 * (long double)cx[j * n + i] * rx[i];
                for (size_t k = j + 1; k < count; ++k)
                    cross += 2 * (long double)cx[j * n + i] * cx[k * n + i];
            }
        }
        CHECK(fabsl(total - diagonal - cross) < 1e-11L);
        if (c.mode == PRISM_BANDS_FIXED)
            CHECK(fabsl(cross) < 1e-11L);
        prism_workspace_destroy(w);
        prism_plan_destroy(p);
    }
}
int main(void) {
    property_test(prism_backend_reference());
    if (prism_backend_fftw())
        property_test(prism_backend_fftw());
    puts("100 seeded linearity, DC, reconstruction and energy cases per backend passed");
    return 0;
}
