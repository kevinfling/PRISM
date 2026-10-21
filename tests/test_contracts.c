#include "check.h"
#include <float.h>
#include <math.h>
#include <prism.h>
#include <stdint.h>
#include <string.h>
typedef struct {
    size_t calls, fail_at, live;
    int frozen;
} allocation_state;
static void *tracked_allocate(void *user, size_t bytes) {
    allocation_state *s = user;
    CHECK(!s->frozen);
    if (s->calls++ == s->fail_at)
        return NULL;
    void *p = malloc(bytes);
    if (p)
        ++s->live;
    return p;
}
static void tracked_free(void *user, void *p) {
    allocation_state *s = user;
    if (p) {
        CHECK(s->live > 0);
        --s->live;
        free(p);
    }
}
static int failed_create(size_t n, const prism_allocator *a, void **s) {
    (void)n;
    (void)a;
    (void)s;
    return PRISM_EBACKEND;
}
static int bad_forward(void *s, const double *x, double *y) {
    (void)s;
    (void)x;
    (void)y;
    return PRISM_EBACKEND;
}
static int nan_forward(void *s, const double *x, double *y) {
    (void)s;
    (void)x;
    y[0] = NAN;
    return PRISM_OK;
}
static void allocation_test(const prism_backend *backend) {
    size_t edges[] = {8, 4};
    prism_config c = {8, 1, 1, PRISM_BANDS_FIXED, edges, NULL};
    int reached_success = 0;
    for (size_t fail = 0; fail < 16; ++fail) {
        allocation_state state = {0, fail, 0, 0};
        prism_allocator a = {&state, tracked_allocate, tracked_free};
        prism_plan *p = NULL;
        prism_workspace *w = NULL;
        int status = prism_plan_create(&c, backend, &a, &p);
        if (status == PRISM_OK)
            status = prism_workspace_create(p, &w);
        if (status == PRISM_OK) {
            double x[8] = {1, 2, 3, 4, 5, 6, 7, 8}, comp[8], res[8], back[8];
            prism_result r = {comp, 8, res, 8, 0, 0, 0, 0};
            state.frozen = 1;
            CHECK(prism_decompose(p, w, x, &r) == PRISM_OK);
            CHECK(prism_reconstruct(&r, back, 8) == PRISM_OK);
            state.frozen = 0;
            reached_success = 1;
        } else
            CHECK(status == PRISM_ENOMEM);
        prism_workspace_destroy(w);
        prism_plan_destroy(p);
        CHECK(state.live == 0);
    }
    CHECK(reached_success);
}
int main(void) {
    size_t edges[] = {8, 4};
    double times[] = {1};
    prism_config c = {8, 1, 1, PRISM_BANDS_FIXED, edges, times};
    prism_plan *p = NULL, *other = NULL;
    prism_workspace *w = NULL;
    CHECK(prism_plan_create(NULL, NULL, NULL, &p) == PRISM_EINVAL);
    CHECK(prism_plan_create(&c, NULL, NULL, NULL) == PRISM_EINVAL);
    c.n = 1;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    c.n = 8;
    c.component_count = SIZE_MAX;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EOVERFLOW);
    c.component_count = 1;
    c.n = SIZE_MAX;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EOVERFLOW);
    c.n = 8;
    c.sample_spacing = NAN;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    c.sample_spacing = 1;
    edges[1] = 8;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    edges[1] = 4;
    c.mode = PRISM_BANDS_HEAT_TELESCOPE;
    times[0] = -1;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    times[0] = INFINITY;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    times[0] = 1;
    c.mode = PRISM_BANDS_HEAT_TELESCOPE;
    c.heat_times = NULL;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    c.heat_times = times;
    times[0] = 0;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    times[0] = 1;
    double repeated[] = {1, 1};
    c.component_count = 2;
    c.heat_times = repeated;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    c.component_count = 1;
    c.heat_times = times;
    c.sample_spacing = DBL_MAX;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_ENUMERIC);
    c.sample_spacing = 1;
    c.mode = PRISM_BANDS_FIXED;
    c.band_edges = NULL;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    c.band_edges = edges;
    edges[1] = 0;
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_EINVAL);
    edges[1] = 4;
    prism_allocator invalid = {NULL, NULL, NULL};
    CHECK(prism_plan_create(&c, NULL, &invalid, &p) == PRISM_EINVAL);
    CHECK(prism_workspace_create(NULL, &w) == PRISM_EINVAL);
    CHECK(prism_plan_create(&c, NULL, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    CHECK(prism_plan_create(&c, NULL, NULL, &other) == PRISM_OK);
    double x[8] = {1, 2, 3, 4, 5, 6, 7, 8}, comp[9], res[9], back[9], weights[16];
    for (size_t i = 0; i < 9; ++i) {
        comp[i] = 1234;
        res[i] = 5678;
    }
    prism_result r = {comp, 7, res, 8, 99, 99, 0, 0};
    CHECK(prism_decompose(p, w, x, &r) == PRISM_ECAPACITY);
    CHECK(comp[0] == 1234 && res[0] == 5678 && r.count == 0);
    r.component_capacity = 8;
    r.residual_capacity = 7;
    CHECK(prism_decompose(p, w, x, &r) == PRISM_ECAPACITY);
    r.residual_capacity = 8;
    CHECK(prism_decompose(other, w, x, &r) == PRISM_EINVAL);
    r.components = x;
    CHECK(prism_decompose(p, w, x, &r) == PRISM_EINVAL);
    r.components = comp;
    r.residual = comp + 1;
    CHECK(prism_decompose(p, w, x, &r) == PRISM_EINVAL);
    r.residual = res;
    x[0] = NAN;
    CHECK(prism_decompose(p, w, x, &r) == PRISM_ENUMERIC);
    x[0] = 1;
    CHECK(prism_plan_weights(p, weights, 15) == PRISM_ECAPACITY);
    CHECK(prism_decompose(p, w, x, &r) == PRISM_OK);
    CHECK(comp[8] == 1234 && res[8] == 5678);
    CHECK(prism_reconstruct(&r, comp, 8) == PRISM_EINVAL);
    CHECK(prism_reconstruct(&r, back, 7) == PRISM_ECAPACITY);
    CHECK(prism_reconstruct(&r, back, 8) == PRISM_OK);
    r.n = SIZE_MAX;
    CHECK(prism_reconstruct(&r, back, 8) == PRISM_EOVERFLOW);
    r.n = 8;
    comp[0] = INFINITY;
    CHECK(prism_reconstruct(&r, back, 8) == PRISM_ENUMERIC);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    prism_plan_destroy(other);
    prism_backend broken = *prism_backend_reference();
    broken.create = failed_create;
    CHECK(prism_plan_create(&c, &broken, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_EBACKEND && !w);
    prism_plan_destroy(p);
    broken = *prism_backend_reference();
    broken.forward = bad_forward;
    CHECK(prism_plan_create(&c, &broken, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    CHECK(prism_decompose(p, w, x, &r) == PRISM_EBACKEND && r.count == 0);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    broken.forward = nan_forward;
    CHECK(prism_plan_create(&c, &broken, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    CHECK(prism_decompose(p, w, x, &r) == PRISM_ENUMERIC);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    broken = *prism_backend_reference();
    broken.inverse = bad_forward;
    CHECK(prism_plan_create(&c, &broken, NULL, &p) == PRISM_OK);
    CHECK(prism_workspace_create(p, &w) == PRISM_OK);
    CHECK(prism_decompose(p, w, x, &r) == PRISM_EBACKEND);
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    if (prism_backend_fftw()) {
        CHECK(prism_plan_create(&c, prism_backend_fftw(), NULL, &p) == PRISM_OK);
        CHECK(prism_workspace_create(p, &w) == PRISM_OK);
        for (size_t i = 0; i < 8; ++i)
            x[i] = DBL_MAX;
        CHECK(prism_decompose(p, w, x, &r) == PRISM_ENUMERIC);
        prism_workspace_destroy(w);
        prism_plan_destroy(p);
    }
    allocation_test(prism_backend_reference());
    if (prism_backend_fftw())
        allocation_test(prism_backend_fftw());
    prism_workspace_destroy(NULL);
    prism_plan_destroy(NULL);
    CHECK(strcmp(prism_status_string(PRISM_OK), "success") == 0);
    puts("capacity, ownership, allocation and backend failure contracts passed");
    return 0;
}
