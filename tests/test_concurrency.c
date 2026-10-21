#include "check.h"
#include <math.h>
#include <prism.h>
#include <pthread.h>
typedef struct {
    const prism_plan *p;
    prism_workspace *w;
    double input[257];
    int failed;
} worker;
static void *run(void *arg) {
    worker *a = arg;
    double comp[514], res[257];
    prism_result r = {comp, 514, res, 257, 0, 0, 0, 0};
    for (size_t i = 0; i < 20; ++i)
        if (prism_decompose(a->p, a->w, a->input, &r) != PRISM_OK ||
            r.reconstruction_relative_error > 1e-10)
            a->failed = 1;
    return NULL;
}
int main(void) {
    const prism_backend *backends[] = {prism_backend_reference(), prism_backend_fftw()};
    for (size_t b = 0; b < 2; ++b) {
        if (!backends[b])
            continue;
        double times[] = {1, 8};
        prism_config c = {257, 2, 1, PRISM_BANDS_HEAT_TELESCOPE, NULL, times};
        prism_plan *p = NULL;
        CHECK(prism_plan_create(&c, backends[b], NULL, &p) == PRISM_OK);
        worker workers[2];
        pthread_t threads[2];
        for (size_t i = 0; i < 2; ++i) {
            workers[i].p = p;
            workers[i].failed = 0;
            CHECK(prism_workspace_create(p, &workers[i].w) == PRISM_OK);
            for (size_t j = 0; j < 257; ++j)
                workers[i].input[j] = sin((double)j * (0.1 + (double)i)) + 1 + (double)i;
        }
        for (size_t i = 0; i < 2; ++i)
            CHECK(pthread_create(&threads[i], NULL, run, &workers[i]) == 0);
        for (size_t i = 0; i < 2; ++i)
            CHECK(pthread_join(threads[i], NULL) == 0);
        for (size_t i = 0; i < 2; ++i) {
            CHECK(!workers[i].failed);
            prism_workspace_destroy(workers[i].w);
        }
        prism_plan_destroy(p);
    }
    puts("concurrent executions with shared plan and separate workspaces passed");
    return 0;
}
