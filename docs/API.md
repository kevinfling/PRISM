# API and ownership

`include/prism.h` is the complete interface. `prism/backend.h`,
`prism/partition.h`, and `prism/heat.h` are convenience includes.

## Lifetimes

`prism_plan_create` validates configuration and builds owned weights. The
configuration and its edge/time arrays may be released after creation. A plan
must outlive its workspaces. Destroy workspaces before destroying the plan.
Failed constructors set the output handle to NULL. NULL destruction is allowed.

An allocator supplies allocation/free callbacks and a borrowed user pointer.
It must return memory aligned for ordinary C types and accept every allocation
made through that allocator. PRISM copies the callbacks. Their code and user
state must remain valid through destruction. FFTW manages its own planner
allocations internally; PRISM-owned memory uses the supplied allocator.

A backend descriptor is borrowed and must remain immutable and alive through
plan destruction. Its create callback produces workspace-local state using the
supplied allocator. That allocator descriptor pointer is valid only during the
create callback; copy it if retaining callbacks or user state. Successful creation
must return non-NULL state. Any partial state returned on failure must be destroyable.
Transforms implement orthonormal DCT-II and transpose synthesis. Execute
callbacks take disjoint arrays, preserve input, write all n outputs, and do not
allocate. PRISM validates finite transform output. Backend error codes propagate.

Plans are read-only after construction. Separate workspaces permit simultaneous
execution. One workspace cannot be shared across overlapping calls. FFTW
workspace creation and destruction use FFTW planner operations and must be
serialized by the caller. The built-in backends support concurrent execution on
separate workspaces. Custom backends must isolate mutable execution state per
workspace and provide the same concurrency contract; the library cannot verify
arbitrary callback behavior.

## Partitions

For J components, fixed edges have J+1 entries:
`n=b0 > b1 > ... > bJ >= 1`. Component j includes indices
`[b(j+1),bj)`. Residual includes `[0,bJ)`. Bands are orthogonal in exact arithmetic.

Heat configuration has J strictly increasing positive finite times.
`Sj[k]=exp(-tj*(pi*k/L)^2)` with S0=1 and L=n*dx. Component j uses S(j)-S(j+1),
and residual uses SJ. The differences are computed with expm1 for small diffusion
increments. Heat bands overlap, so energy accounting includes cross terms.
Samples are cell-centered; times have length-squared units. This is the continuum
Neumann spectral convention, not a finite-difference stencil eigenvalue.

## Execution and results

Input has n finite doubles. Components have J*n doubles, in high-to-low order;
residual has n doubles. Capacity fields count doubles. Each actual used buffer
region must be disjoint, including partial overlap. No in-place operation is
supported. Unused capacity is untouched. Buffers are caller-owned and need no
FFTW-specific alignment. All size products are checked before allocation.

`prism_decompose` sets count/n and diagnostics to zero on entry. Invalid arguments,
insufficient capacities, overlap, and nonfinite input are rejected before writing
signal output buffers. Backend/numerical failures may leave partial contents;
only consume results after PRISM_OK. Successful count/n describe valid outputs.

Diagnostics use long-double accumulation and hypot norms. Relative error is
`norm(x-reconstruction)/norm(x)` for nonzero input; zero input reports absolute
error. Max absolute error is recorded separately. A finite input whose transform
or synthesized result exceeds representable double range returns PRISM_ENUMERIC.
The 1e-10 target applies to tested well-scaled inputs, not arbitrary extremes.
PRISM_OK does not enforce that target or the partition-defect target. Callers may
inspect the diagnostics and choose an application-specific acceptance threshold.
Diagnostics accumulate stored outputs in long double, whose precision is
platform-dependent; they are not certified error bounds. The reconstructed
double buffer incurs an additional final rounding.

`prism_reconstruct` sums components and residual into a separate caller buffer.
`prism_plan_weights` copies J+1 contiguous rows, with residual weights last.
`prism_status_string` describes status codes. Custom callbacks should return only
PRISM status values. PRISM's explicit allocations occur during plan/workspace
construction only.
Custom backend execution callbacks must not allocate. The allocator test observes
allocations through the provided callbacks; it does not interpose on every
possible allocation inside external libraries.

## Release and error reference

`PRISM_VERSION` and its MAJOR/MINOR/PATCH macros identify the header release.
The initial API is pre-1.0. CMake compatibility and shared-library versioning are
scoped to the 0.1 series; later minor versions may change API or ABI.

| Status | Checked condition |
| --- | --- |
| `PRISM_OK` | Operation completed; decomposition outputs and diagnostics are finite |
| `PRISM_EINVAL` | Invalid arguments, unsupported mode, incompatible workspace, or overlapping signal buffers |
| `PRISM_ECAPACITY` | A supplied output capacity is too small |
| `PRISM_ENOMEM` | A PRISM-owned allocation failed |
| `PRISM_EBACKEND` | Backend construction/execution failed or returned invalid state |
| `PRISM_ENUMERIC` | Nonfinite input, nonfinite transform/output/diagnostic, or a checked numerical range failure |
| `PRISM_EOVERFLOW` | A checked size product or backend dimension exceeds its supported range |

A backend may propagate any PRISM error code. The built-in FFTW workspace limits
n to INT_MAX because its one-dimensional planner accepts an int. Ordinary C
buffer validity remains a caller responsibility: capacities describe allocation
sizes truthfully, pointers designate live appropriately aligned double arrays,
and config arrays contain the declared number of entries. All dimensions, strides,
and lifetimes must agree. The result descriptor itself must remain disjoint from
signal buffers, and custom callbacks must satisfy their contracts.

`prism_reconstruct` validates signal-buffer disjointness and capacities. Numerical
failure can leave a partially written destination, which must not be consumed.
`prism_plan_weights` needs capacity (J+1)*n doubles and copies the residual row last.
The convenience headers add no separate modules or algorithms.
