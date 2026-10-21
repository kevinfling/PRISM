/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_NOISE_H
#define PRISM_NOISE_H
#include <prism.h>
#include <prism/rng.h>
#ifdef __cplusplus
extern "C" {
#endif

typedef struct prism_noise_plan prism_noise_plan;
typedef struct prism_noise_workspace prism_noise_workspace;
typedef enum {
    PRISM_NOISE_CHOLESKY,
    PRISM_NOISE_CIRCULANT_REFERENCE,
    PRISM_NOISE_CIRCULANT_FFTW
} prism_noise_method;
typedef struct {
    size_t n;
    prism_noise_method method;
    /* Circulant methods only: [0,1). Zero rejects every negative computed
     * eigenvalue. Positive values explicitly permit rounding small negatives
     * to zero relative to the largest absolute computed eigenvalue. Such a
     * correction changes the target covariance and is reported below. */
    double embedding_relative_tolerance;
} prism_noise_config;
typedef struct {
    prism_noise_method method;
    size_t n, embedding_size, rounded_eigenvalue_count;
    double minimum_embedding_eigenvalue, maximum_embedding_eigenvalue;
    double maximum_embedding_imaginary_part, maximum_eigenvalue_adjustment;
    double minimum_cholesky_pivot;
} prism_noise_diagnostics;

int prism_noise_method_available(prism_noise_method);
/* Unit-spaced fGn increment covariance at lags 0..count-1, scaled by variance.
 * 0 < hurst < 1; variance >= 0; count >= 1. This is fGn, not fBm levels.
 * Sampling interval scaling, if needed, belongs in the explicit variance. */
int prism_fgn_covariance(size_t count, double hurst, double variance, double *out, size_t capacity);
/* Toeplitz target K[i,j] = covariance_lags[abs(i-j)]. Arrays are copied.
 * Cholesky uses the dense Toeplitz matrix. Circulant methods use a fixed mirror
 * embedding of length 2*(n-1), or 1 for n=1. A PSD Toeplitz target need not have
 * a PSD embedding; rejection does not establish that the target is invalid.
 * Diagnostics are mandatory, initialized even on failure, and method-specific
 * fields not applicable to a method are zero. Failed creation sets *out=NULL.
 * Descriptor objects (config, allocator, handles, diagnostics) and lag arrays
 * must be mutually disjoint. Serialize FFTW plan/workspace creation and destruction with other
 * planners. */
int prism_noise_plan_create(const prism_noise_config *, const double *covariance_lags,
                            size_t covariance_capacity, const prism_allocator *,
                            prism_noise_plan **out, prism_noise_diagnostics *);
void prism_noise_plan_destroy(prism_noise_plan *);
int prism_noise_workspace_create(const prism_noise_plan *, prism_noise_workspace **out);
void prism_noise_workspace_destroy(prism_noise_workspace *);
/* No execute-time allocation. Separate workspaces and RNGs allow concurrent
 * sampling. Output must not overlap RNG state. Validate all capacities before
 * consuming RNG state. Numerical failure may consume RNG state, but leaves the
 * caller's output unchanged. Plans must outlive workspaces. */
int prism_noise_sample(const prism_noise_plan *, prism_noise_workspace *, prism_rng *, double *out,
                       size_t capacity);
/* Number of white innovations required: n for dense, embedding length for
 * circulant methods. NULL plan returns zero. No RNG or allocation is involved. */
size_t prism_noise_innovation_count(const prism_noise_plan *);
/* Apply the covariance square root to caller-provided finite innovations.
 * This operation determines neither their distribution nor independence.
 * Input/output must be disjoint. Failures leave output unchanged. */
int prism_noise_transform(const prism_noise_plan *, prism_noise_workspace *,
                          const double *innovations, size_t innovation_capacity, double *out,
                          size_t capacity);
/* Custom-source variant of prism_noise_sample. Same transactional output and
 * workspace contract. Descriptor, state, and output must be disjoint. Errors
 * propagate from the callback; nonfinite innovations return PRISM_ENUMERIC. */
int prism_noise_sample_from(const prism_noise_plan *, prism_noise_workspace *,
                            const prism_normal_source *, double *out, size_t capacity);
/* Copy the requested target lags, not the covariance after optional rounding. */
int prism_noise_target_covariance(const prism_noise_plan *, double *out, size_t capacity);
#ifdef __cplusplus
}
#endif
#endif
