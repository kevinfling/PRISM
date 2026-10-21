#include "check.h"
#include <float.h>
#include <math.h>
#include <prism.h>
#include <stdint.h>
#include <string.h>
static double largest_relative_error = 0;
static const double pi = 3.14159265358979323846;
static void *allocate(void *user, size_t n) {
    (void)user;
    return malloc(n);
}
static void deallocate(void *user, void *p) {
    (void)user;
    free(p);
}
/* Independent direct cosine matrix oracle; no library implementation reused. */
static void oracle(const double *x, double *y, size_t n, int inverse) {
    for (size_t a = 0; a < n; ++a) {
        long double sum = 0;
        for (size_t b = 0; b < n; ++b) {
            size_t k = inverse ? b : a, i = inverse ? a : b;
            double basis = sqrt((k ? 2.0 : 1.0) / (double)n) *
                           cos(pi * ((double)i + 0.5) * (double)k / (double)n);
            sum += (long double)basis * (long double)x[b];
        }
        y[a] = (double)sum;
    }
}
static void transform_test(const prism_backend *b) {
    const size_t n = 31;
    double x[31], y[31], expected[31], back[31];
    prism_allocator a = {NULL, allocate, deallocate};
    void *state = NULL;
    CHECK(b->create(n, &a, &state) == PRISM_OK);
    for (int kind = 0; kind < 4; ++kind) {
        for (size_t i = 0; i < n; ++i)
            x[i] = kind == 0   ? 1.0
                   : kind == 1 ? (i == n / 2 ? 1.0 : 0.0)
                   : kind == 2 ? cos(pi * ((double)i + 0.5) * 7.0 / (double)n)
                               : sin((double)i * 0.91) + 0.03 * (double)i;
        oracle(x, expected, n, 0);
        CHECK(b->forward(state, x, y) == PRISM_OK);
        for (size_t i = 0; i < n; ++i)
            CHECK(fabs(y[i] - expected[i]) < 1e-12);
        oracle(y, expected, n, 1);
        CHECK(b->inverse(state, y, back) == PRISM_OK);
        for (size_t i = 0; i < n; ++i) {
            CHECK(fabs(back[i] - expected[i]) < 1e-12);
            CHECK(fabs(back[i] - x[i]) < 1e-12);
        }
    }
    /* All DCT basis vectors, including DC and the highest frequency. */
    for (size_t k = 0; k < n; ++k) {
        memset(x, 0, sizeof(x));
        x[k] = 1;
        oracle(x, expected, n, 1);
        CHECK(b->inverse(state, x, y) == PRISM_OK);
        for (size_t i = 0; i < n; ++i)
            CHECK(fabs(y[i] - expected[i]) < 1e-12);
    }
    b->destroy(state);
}
static void decompose_test(const prism_backend *b, size_t n, prism_mode mode) {
    size_t edges[3] = {n, n / 2, n / 8};
    double times[2] = {0.5, 8.0};
    prism_config c = {n, 2, 0.25, mode, edges, times};
    prism_plan *p = NULL;
    prism_workspace *w = NULL, *other = NULL;
    CHECK(prism_plan_create(&c, b, NULL, &p) == PRISM_OK);
    /* Configuration is consumed at plan creation. */
    edges[1] = 0;
    times[0] = -1;
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    CHECK(prism_workspace_create(p, &other) == PRISM_OK);
    double *x = malloc(n * sizeof(double)), *components = malloc(2 * n * sizeof(double));
    double *residual = malloc(n * sizeof(double)), *reconstructed = malloc(n * sizeof(double));
    double *weights = malloc(3 * n * sizeof(double));
    CHECK(x && components && residual && reconstructed && weights);
    prism_result r = {components, 2 * n, residual, n, 0, 0, 0, 0};
    CHECK(prism_plan_weights(p, weights, 3 * n) == PRISM_OK);
    for (size_t k = 0; k < n; ++k) {
        double sum = weights[k] + weights[n + k] + weights[2 * n + k];
        CHECK(fabs(sum - 1) < 1e-14);
        for (size_t j = 0; j < 3; ++j)
            CHECK(weights[j * n + k] >= 0 && weights[j * n + k] <= 1);
        if (mode == PRISM_BANDS_HEAT_TELESCOPE) {
            double omega = pi * (double)k / ((double)n * 0.25);
            double s1 = exp(-0.5 * omega * omega), s2 = exp(-8.0 * omega * omega);
            CHECK(fabs(weights[k] - (1 - s1)) < 1e-14);
            CHECK(fabs(weights[n + k] - (s1 - s2)) < 1e-14);
            CHECK(fabs(weights[2 * n + k] - s2) < 1e-14);
        } else {
            CHECK(weights[k] == (k >= n / 2 ? 1.0 : 0.0));
            CHECK(weights[n + k] == (k >= n / 8 && k < n / 2 ? 1.0 : 0.0));
        }
    }
    for (int trial = 0; trial < (n > 1024 ? 1 : 7); ++trial) {
        double scale = trial == 3 ? 1e-150 : trial == 4 ? 1e150 : 1;
        for (size_t i = 0; i < n; ++i) {
            double signal =
                trial == 0 ? sin(0.19 * (double)i) + 0.7 * cos(0.71 * (double)i) + 0.01 * (double)i
                : trial == 1 ? 0
                : trial == 2 ? 3.0
                : trial == 5 ? (i == n - 1 ? 1.0 : 0.0)
                : trial == 6 ? cos(pi * ((double)i + 0.5) * (double)(n - 1) / (double)n)
                             : sin(0.37 * (double)i);
            x[i] = scale * signal;
        }
        CHECK(prism_decompose(p, trial % 2 ? w : other, x, &r) == PRISM_OK);
        CHECK(r.count == 2 && r.n == n);
        CHECK(r.reconstruction_relative_error <= 1e-10);
        if (r.reconstruction_relative_error > largest_relative_error)
            largest_relative_error = r.reconstruction_relative_error;
        CHECK(prism_reconstruct(&r, reconstructed, n) == PRISM_OK);
        for (size_t i = 0; i < n; ++i)
            CHECK(fabs((reconstructed[i] - x[i]) / scale) < 1e-9);
        /* Oracle checks components, so reconstruction alone cannot hide bad bands. */
        if (n == 128 && trial == 0) {
            double spectrum[128], filtered[128], expected[128];
            oracle(x, spectrum, n, 0);
            for (size_t j = 0; j < 3; ++j) {
                for (size_t k = 0; k < n; ++k)
                    filtered[k] = spectrum[k] * weights[j * n + k];
                oracle(filtered, expected, n, 1);
                const double *actual = j == 2 ? residual : components + j * n;
                for (size_t i = 0; i < n; ++i)
                    CHECK(fabs(expected[i] - actual[i]) < 1e-11);
            }
        }
    }
    free(weights);
    free(reconstructed);
    free(residual);
    free(components);
    free(x);
    prism_workspace_destroy(other);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
}
int main(int argc, char **argv) {
    const prism_backend *b =
        argc > 1 && strcmp(argv[1], "fftw") == 0 ? prism_backend_fftw() : prism_backend_reference();
    CHECK(b);
    transform_test(b);
    const size_t sizes[] = {128, 257, 511, 512, 513, 1024, 4096, 16384};
    for (size_t i = 0; i < sizeof(sizes) / sizeof(sizes[0]); ++i) {
        decompose_test(b, sizes[i], PRISM_BANDS_FIXED);
        decompose_test(b, sizes[i], PRISM_BANDS_HEAT_TELESCOPE);
        printf("%s n=%zu passed\n", b->name, sizes[i]);
        fflush(stdout);
    }
    printf("%s largest relative reconstruction error: %.9g\n", b->name, largest_relative_error);
    return 0;
}
