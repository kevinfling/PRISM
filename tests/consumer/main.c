#include <math.h>
#include <prism.h>
#include <prism/noise.h>
#include <prism/sigma.h>
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
    if (!status) {
        prism_rng rng;
        double lags[4], sample[4];
        prism_noise_plan *np = 0;
        prism_noise_workspace *nw = 0;
        prism_noise_diagnostics diag;
        prism_noise_config nc = {4, PRISM_NOISE_CHOLESKY, 0};
        status = prism_fgn_covariance(4, 0.7, 1, lags, 4);
        if (!status)
            status = prism_noise_plan_create(&nc, lags, 4, 0, &np, &diag);
        if (!status)
            status = prism_noise_workspace_create(np, &nw);
        if (!status)
            status = prism_rng_seed(&rng, 42, 54);
        if (!status) {
            prism_normal_source source = prism_rng_source(&rng);
            status = prism_noise_sample_from(np, nw, &source, sample, 4);
        }
        if (!status && prism_sigma_available()) {
            prism_sigma_rng srng;
            status = prism_sigma_seed(&srng, 42, 54);
            prism_normal_source source = prism_sigma_source(&srng);
            if (!status)
                status = prism_noise_sample_from(np, nw, &source, sample, 4);
        }
        prism_noise_workspace_destroy(nw);
        prism_noise_plan_destroy(np);
    }
    return status;
}
