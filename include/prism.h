/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Kevin Fling */
#ifndef PRISM_H
#define PRISM_H
#include <stddef.h>
#define PRISM_VERSION_MAJOR 0
#define PRISM_VERSION_MINOR 1
#define PRISM_VERSION_PATCH 0
#define PRISM_VERSION "0.1.0"
#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    PRISM_OK = 0,
    PRISM_EINVAL,
    PRISM_ECAPACITY,
    PRISM_ENOMEM,
    PRISM_EBACKEND,
    PRISM_ENUMERIC,
    PRISM_EOVERFLOW
} prism_status;
typedef enum { PRISM_BANDS_FIXED, PRISM_BANDS_HEAT_TELESCOPE } prism_mode;
typedef struct prism_plan prism_plan;
typedef struct prism_workspace prism_workspace;
typedef struct {
    void *user;
    void *(*allocate)(void *user, size_t bytes);
    void (*deallocate)(void *user, void *ptr);
} prism_allocator;

/* Backends implement an orthonormal DCT-II and its transpose. create/destroy
 * own workspace-local state; execution must not allocate or mutate inputs.
 * Every output is disjoint from its input. A backend descriptor must remain
 * immutable and outlive its plans. create must copy the temporary allocator
 * descriptor if retaining it, and return non-NULL state on success. Callbacks
 * return a prism_status; numerical output is checked by prism_decompose.
 * Mutable execution state must be isolated per workspace. */
typedef struct {
    const char *name;
    int (*create)(size_t n, const prism_allocator *allocator, void **state);
    int (*forward)(void *state, const double *input, double *output);
    int (*inverse)(void *state, const double *input, double *output);
    void (*destroy)(void *state);
} prism_backend;
const prism_backend *prism_backend_reference(void);
/* Returns NULL when PRISM was built without FFTW. */
const prism_backend *prism_backend_fftw(void);
const char *prism_fftw_version(void);

typedef struct {
    size_t n, component_count;
    double sample_spacing;
    prism_mode mode;
    const size_t *band_edges; /* fixed: count+1 descending edges, first=n */
    const double *heat_times; /* heat: count increasing positive times */
} prism_config;
typedef struct {
    double *components;        /* row-major [count][n] */
    size_t component_capacity; /* doubles */
    double *residual;
    size_t residual_capacity; /* doubles */
    size_t count, n;          /* written on success */
    double reconstruction_relative_error, reconstruction_max_abs_error;
} prism_result;

/* NULL backend selects FFTW if available, otherwise reference. NULL allocator
 * selects malloc/free. Configuration arrays are copied into partition weights.
 * Allocator user state must outlive the plan and all its workspaces. */
int prism_plan_create(const prism_config *, const prism_backend *, const prism_allocator *,
                      prism_plan **);
void prism_plan_destroy(prism_plan *);
int prism_workspace_create(const prism_plan *, prism_workspace **);
void prism_workspace_destroy(prism_workspace *);
/* Plan must outlive its workspaces. Separate workspaces permit concurrent
 * execution. Serialize FFTW workspace creation/destruction across threads.
 * Input must be finite; input and output buffers must be disjoint. Capacities
 * are checked before writes. Numerical/backend failures may leave partial buffer
 * contents; result count/n are zero on failure. PRISM's execute paths perform no
 * explicit allocations. PRISM_OK does not enforce a reconstruction tolerance;
 * callers can apply their own acceptance criteria to the reported diagnostics. */
int prism_decompose(const prism_plan *, prism_workspace *, const double *, prism_result *);
int prism_reconstruct(const prism_result *, double *output, size_t capacity);
/* Copies count+1 rows of spectral weights (residual last), for diagnostics. */
int prism_plan_weights(const prism_plan *, double *output, size_t capacity);
const char *prism_status_string(int status);
#ifdef __cplusplus
}
#endif
#endif
