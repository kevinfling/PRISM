#include <math.h>
#include <prism.h>
int main(void) {
    size_t edges[] = {4, 2};
    prism_config c = {4, 1, 1, PRISM_BANDS_FIXED, edges, 0};
    double x[] = {1, 2, 3, 4}, components[4], residual[4];
    prism_result r = {components, 4, residual, 4, 0, 0, 0, 0};
    prism_plan *p = 0;
    prism_workspace *w = 0;
    int status = prism_plan_create(&c, 0, 0, &p);
    if (!status)
        status = prism_workspace_create(p, &w);
    if (!status)
        status = prism_decompose(p, w, x, &r);
    if (!status && r.reconstruction_relative_error > 1e-10)
        status = 1;
    prism_workspace_destroy(w);
    prism_plan_destroy(p);
    return status;
}
