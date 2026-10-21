# Covariance sampling and the fixed-linear experiment

This v0.2.0 development module supplies explicit seeded sampling and a research
experiment. Deterministic decomposition remains the default. Sampling does not
estimate a covariance, select a partition, calibrate uncertainty, or demonstrate
application benefit. Include [prism/noise.h](../include/prism/noise.h) explicitly;
[the complete example](../examples/noise.c) uses the dependency-free dense method.

## Covariance definition and scale

For a real stationary target, supply n finite lags and define

```math
K_{ij}=\gamma(|i-j|),\qquad 0\le i,j\lt n.
```

A valid Gaussian covariance must be positive semidefinite. Nonnegative lag-zero
variance and $`|\gamma(k)|\le\gamma(0)`$ are necessary but insufficient. The selected
factorization performs additional numerical checks; failure does not constitute
a mathematical proof that a near-singular target is invalid.

`prism_fgn_covariance` supplies the covariance of unit-spaced increments of
fractional Brownian motion with Hurst parameter $`0\lt H\lt1`$ and explicit
increment variance $`v\ge0`$:

```math
\gamma_H(0)=v,\qquad
\gamma_H(k)=\frac v2\left[(k+1)^{2H}-2k^{2H}+(k-1)^{2H}\right],\quad k\ge1.
```

These are stationary **fGn increments**, rather than fBm levels. At H=1/2 the
implementation sets all nonzero lags exactly to zero. With an fBm scale a and
sampling interval h, use $`v=a^2h^{2H}`$ explicitly; PRISM does not infer sampling
units. See Dieker and Mandjes, [definitions and simulation methods, §2](https://www.columbia.edu/~ad3217/publications/specsim.pdf).

Direct subtraction of nearly equal powers is unstable for large k or H near
1/2. Let p=2H. At k=1 the normalized covariance is
$`\mathop{\mathrm{expm1}}((p-1)\log2)`$. At k≥2 the code evaluates the even
binomial expansion

```math
\frac{\gamma_H(k)}v=
\frac{p(p-1)}2 k^{p-2}\sum_{j=0}^{\infty}a_j,
\qquad a_0=1,\qquad
\frac{a_{j+1}}{a_j}=
\frac{(p-2j-2)(p-2j-3)}{(2j+3)(2j+4)k^2}.
```

This follows by expanding both $`(1+1/k)^p`$ and $`(1-1/k)^p`$; odd powers cancel.
All normalized series terms are positive for 0<p<2, with successive ratios
below 1/4 when k≥2. The implementation accumulates in long double until its
relative epsilon criterion or 127 additional terms. This avoids that specific
cancellation; it is not a certified relative-error bound, and long-double
precision is platform-dependent. Variance zero is handled without division.

## Dense reference factorization

`PRISM_NOISE_CHOLESKY` forms a lower triangular factor with the recurrence

```math
L_{ii}=\sqrt{K_{ii}-\sum_{k\lt i}L_{ik}^2},\qquad
L_{ij}=\frac{K_{ij}-\sum_{k\lt j}L_{ik}L_{jk}}{L_{jj}},\quad j\lt i.
```

The inner sums use long double and stored factors use double. A negative computed
pivot is rejected; a zero pivot requires a zero off-diagonal remainder. There
is no implicit jitter, eigenspectrum projection, or pivot tolerance. Zero and
some singular targets are supported, including the tested constant unit
covariance, but arbitrary numerically singular targets need not be accepted.
`minimum_cholesky_pivot` is the smallest computed pre-square-root diagonal
remainder encountered. For independent standard normal z, the exact-arithmetic
construction $`y=Lz`$ has covariance $`LL^{\mathsf T}=K`$.

Plan setup uses O(n³) operations and n² doubles for the factor; sampling uses
O(n²) operations and workspace arrays of length n. This is a small-input
reference, rather than a large-scale default.

## Circulant embedding, normalization, and rejection

Both `PRISM_NOISE_CIRCULANT_REFERENCE` and `PRISM_NOISE_CIRCULANT_FFTW` use
m=2(n−1), with m=1 for n=1. The mirror vector is

```math
c=(\gamma(0),\ldots,\gamma(n-1),\gamma(n-2),\ldots,\gamma(1)).
```

For n=1 use only the lag-zero entry; n=2 has no appended mirror tail.
Let $`F_{kj}=\exp(-2\pi\mathrm{i}kj/m)`$, $`\lambda=Fc`$, and
$`C=F^*\mathop{\mathrm{diag}}(\lambda)F/m`$. This real symmetric circulant
matrix has K as its upper-left n×n block. If all eigenvalues are nonnegative,

```math
B=\frac1m F^*\mathop{\mathrm{diag}}(\sqrt\lambda)F,
\qquad y=Bz,\qquad z\sim N(0,I_m),
\qquad BB^{\mathsf T}=C.
```

The last identity uses $`FF^*=mI`$. Return the first n entries of y. A real
forward transform, multiplication by square-root eigenvalues, and an inverse
transform divided by m implement this construction. DC and even-m Nyquist
coefficients have multiplicity one; interior complex-conjugate pairs have
multiplicity two. There is no extra factor of two in the spectral multiplier.
FFTW's forward/backward transforms are unnormalized ([FFTW definition](https://www.fftw.org/doc/The-1d-Real_002ddata-DFT.html)).
The dependency-free implementation uses direct long-double trigonometric sums,
not a DCT approximation of the Toeplitz matrix. Its O(m²) transforms are intended
for reference checks. FFTW provides the fast transform path with O(m log m)
algorithmic cost, O(m) PRISM storage, and planner setup outside execution.

A PSD target can have an indefinite fixed mirror embedding. For example,
$`\gamma=(1.05,0.9,0.62)`$ gives a positive definite 3×3 Toeplitz target: the
covariance from $`\cos(k\theta)`$, with $`\cos\theta=0.9`$, plus 0.05I.
Its 4×4 embedding has eigenvalues 3.47, 0.43, −0.13, 0.43 and is rejected.
There is no automatic embedding enlargement or fallback to another method.

The mandatory diagnostics include embedding size, minimum/maximum eigenvalues,
and maximum absolute imaginary part of the computed spectrum. Imaginary parts
should vanish analytically; they are reported as a roundoff observation and are
not a certified error bound. The real parts define the numerical factors.

`embedding_relative_tolerance=0` rejects every negative computed eigenvalue.
A caller can explicitly choose $`0\lt\tau\lt1`$ to permit setting a negative
computed eigenvalue to zero only when

```math
-\lambda_k\le\tau\max_j|\lambda_j|.
```

Diagnostics record the full-spectrum multiplicity count of permitted corrections
and the largest correction. With correction $`\delta\lambda_k`$, the effective
circulant covariance is changed by
$`\Delta C=F^*\mathop{\mathrm{diag}}(\delta\lambda)F/m`$ and
$`\|\Delta C\|_2=\max_k|\delta\lambda_k|`$ in exact arithmetic. Cropping cannot
increase that norm. This does not include all floating-point errors.
`prism_noise_target_covariance` copies the **requested** lags, not this modified
covariance. A corrected plan must not be described as sampling the exact target.
Dense Cholesky requires this embedding-only tolerance to be zero.

“Exact covariance” describes these constructions under ideal independent
continuous Gaussian inputs and exact arithmetic. Actual samples use finite-state
pseudorandom numbers, discretized uniforms, libm, and floating-point transforms.
The tests check selected sizes and H values; they establish no all-input or
all-platform guarantee.

## RNG, ownership, errors, and concurrency

[prism/rng.h](../include/prism/rng.h) exposes caller-owned state. Seed it before
use with a 64-bit seed and a stream ID fitting in 63 bits. Integer generation uses
PCG-XSH-RR 64/32; the published seed-42/stream-54 vector is a regression test.
See [the PCG author documentation](https://www.pcg-random.org/using-pcg-c-basic.html).
There is no clock seeding or process-global RNG.

Uniform generation concatenates two successive 32-bit words, takes the high 53
bits, and multiplies by 2⁻⁵³ to produce [0,1) doubles. Box–Muller uses
$`R=\sqrt{-2\log(1-U)}`$ and $`\Theta=2\pi V`$, returning $`R\cos\Theta`$ and
caching $`R\sin\Theta`$. Thus the logarithm receives a value in (0,1]. Chunked
normal draws preserve the sequence, including the cached spare. The finite grid
also bounds the generated radius; these samples are not continuous Gaussian
random variables in a literal mathematical sense. This RNG is for simulation,
not cryptographic use.

The integer sequence is specified exactly; normal values are not promised to be
bitwise identical across libm implementations or platforms. Distinct stream IDs
provide reproducible indexing, rather than a proof of independent random
variables. Assign IDs to logical replicates, not worker IDs. Shared immutable
plans with separate workspaces and RNGs permit concurrent execution. Tests now
compare serial scheduling with 1,2,4,7 threads using stream ID equal to logical
sample ID. Reusing a continuing RNG still requires a fixed draw order within
that source; source indexing does not fix parallel reduction order.

Plans copy lag arrays and allocator callbacks; allocator state must outlive the
plan and its workspaces. Plans must outlive workspaces. NULL destruction is
allowed. Serialize FFTW plan **and** workspace creation/destruction across all
FFTW planner users. FFTW workspace state is private to that workspace.

Constructor descriptors, output handles, diagnostics, and lag arrays must be
mutually disjoint. Capacities count doubles, must truthfully describe live arrays, and are checked
before RNG consumption or output writes. RNG state and output arrays must be
disjoint. Constructors set output handles to NULL before failure; diagnostics
are mandatory and initialized, with inapplicable fields zero. On a failed
factorization they can describe the partial check, rather than a usable plan.
Sampling validates its workspace owner. Numerical sampling failure can consume
RNG state but leaves the caller's output unchanged. No PRISM execution-time
allocation occurs. FFTW planner allocations are outside the caller allocator;
allocation tests do not interpose on arbitrary external-library allocations.

Invalid dimensions, nonfinite lags, tolerance, or RNG parameters return
`PRISM_EINVAL`; short arrays return `PRISM_ECAPACITY`; checked size/int-dimension overflow returns
`PRISM_EOVERFLOW`; missing FFTW returns `PRISM_EBACKEND`. Factorization or embedding
rejection (including a negative variance or a lag exceeding it in magnitude)
returns `PRISM_ECOVARIANCE`, allocation failure `PRISM_ENOMEM`, and
nonfinite numerical results `PRISM_ENUMERIC`. The core decomposition's existing
error and partial-output contract is unchanged.

## Innovation providers and deterministic covariance transforms

[prism/innovation.h](../include/prism/innovation.h) defines `prism_normal_source`,
a borrowed state pointer, its live byte extent, and a checked fill callback.
`prism_rng_source` wraps the built-in PCG state. The optional
[SIGMA adapter](SIGMA.md) supplies another source; applications can provide their
own. The callback must not allocate, must fill the requested count with finite
values on success, and must return PRISM status codes. All reachable state must
remain disjoint from outputs. The descriptor, direct state extent, and output
are checked for overlap before invoking the callback. The normal-fill function
can leave partial output and consume state on callback/numerical failure.

`prism_noise_sample_from` stages those innovations in its workspace before
applying the covariance square root and committing the sample. It preserves
caller output on every failure. The original `prism_noise_sample` convenience
function uses that same path with a PCG source; its random-word consumption and
normal sequence are unchanged. A finite-value check cannot establish that a
custom source is independent, centered, standardized, or Gaussian.

`prism_noise_innovation_count` reports the required white input count: n for
dense Cholesky and m for circulant methods; NULL returns zero.
`prism_noise_transform` accepts that many caller-supplied innovations and applies
the square root without a source, RNG consumption, allocation, or planning.
Innovation input and sample output must be disjoint. Invalid capacities,
nonfinite innovations, or numerical failure leave the output unchanged.

For a general supplied vector z with mean μ and covariance V and implemented
square-root/cropping operator R, exact arithmetic gives

```math
\mathbb E[Rz]=R\mu,\qquad
\mathop{\mathrm{Cov}}(Rz)=RVR^{\mathsf T}.
```

For an uncorrected exact factor, appropriate standardized innovations give
$`RR^{\mathsf T}=K`$; optional eigenvalue correction changes that target as
specified above. Gaussianity is additionally needed to call the resulting
field Gaussian. This boundary keeps innovation generation, covariance
transformation, and later marginal/nonlinear transformations separate. Marginal
transformations and their covariance recalibration are not implemented here.

## Fixed-linear protocol and interpretation

The [research driver](../experiments/ensemble.c) checks the conditional identity
in [README §10](../README.md#10-the-fixed-linear-ensemble-identity-research-only).
It is enabled separately on POSIX systems:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPRISM_BUILD_EXPERIMENTS=ON \
    -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel 3
./build/experiments/prism_ensemble_experiment fftw > ensemble.csv 2> ensemble.txt
./build/experiments/prism_ensemble_experiment reference > ensemble-reference.csv 2> ensemble-reference.txt
# An optional final argument selects the enabled SIGMA adapter instead of PCG:
# ./build/experiments/prism_ensemble_experiment fftw sigma > ensemble-sigma.csv 2> ensemble-sigma.txt
python3 tools/summarize_ensemble.py ensemble.csv
```

The fixed protocol is n=64, dx=1, heat time t=16, q=3, amplitude σ=0.25,
master seed 20261021, and 20 seeded replicates per control. Stream ID is
4×replicate+control. Controls are zero, white, fGn H=0.3, and fGn H=0.7; nonzero
controls have unit marginal variance. Counts are E=1,4,16,64,256, using nested
prefixes within each replicate. Counts therefore have correlated observations;
they must not be treated as independent experimental groups. Parameters are fixed
in source before the archived run, with no application tuning or held-out data.

Writing Q for the README's orthonormal DCT matrix C, define

```math
w_k=(1-e^{-16(\pi k/64)^2})^3,\qquad
T=Q^{\mathsf T}\mathop{\mathrm{diag}}(w)Q,\qquad
x_i=\cos(2\pi(i+\tfrac12)/64)+0.3\cos(23\pi(i+\tfrac12)/64)+0.01i.
```

Each trial computes $`T(x+\sigma\eta_e)`$. The deterministic baseline is Tx;
this is a numerical identity experiment with **no physical ground-truth signal**.
For each replicate, the reported squared error is $`D_E=\|z_E-Tx\|_2^2`$.
Let $`A=TKT^{\mathsf T}`$. Independent ideal Gaussian trials give

```math
\mathbb E D_E=\frac{\sigma^2}E\mathop{\mathrm{tr}}(A),\qquad
\mathop{\mathrm{Var}}(D_E)=\frac{2\sigma^4}{E^2}\mathop{\mathrm{tr}}(A^2).
```

To derive the variance, diagonalize symmetric A with eigenvalues aᵢ. Then D is
a weighted sum of independent squared standard normals with coefficients
$`\sigma^2a_i/E`$; each squared normal has variance two. The driver builds T and
A from a direct cosine matrix independently of the transform implementations.
It checks the paired linearity defect on every replicate and compares the mean
of 20 squared errors with prediction using a six-theoretical-standard-error
smoke threshold. These are fixed-seed checks, not formal significance tests or
multiple-comparison-adjusted confidence statements.

CSV stores all 400 replicate/count observations, theoretical moments, stream
IDs, parameters, and cumulative sampling, transform, and aggregation times.
Setup and one deterministic execution are recorded separately. Independent
oracle/diagnostic costs and printing are outside those component timings;
the rows are not complete end-to-end latency measurements. Metadata names the
backend, noise method, version, FFTW version, and planner policy. The FFTW path
uses circulant sampling; the reference path uses dense Cholesky, so identical
seeds do not imply paired noise realizations across those methods. Their timings
cannot isolate transform-backend performance. Hardware conditions are uncontrolled.

The summary reports mean squared error, empirical standard error across
replicates, and an approximate 95% Student-t interval with 19 degrees of freedom.
Those intervals assume effectively independent replicates; squared-error
observations are not normally distributed, so small-sample coverage is not
certified. They are descriptive simulation uncertainty, not uncertainty in
physical source estimates. For exactly zero errors the reported interval is
zero-width; that is a control, not a calibration claim.

The observed 1/E variance reduction is relative to the noisy single-trial
procedure. Deterministic Tx already has zero perturbation error. No ensemble
benefit over that baseline is asserted. A production ensemble decomposition API,
adaptive stopping, nonlinear partitions, other backend adapters, and auxiliary
operators remain outside this scope.

## Attribution

The PCG state transition, output permutation, and seeding in `src/rng.c` are
adapted from Melissa O'Neill's [minimal C implementation](https://github.com/imneme/pcg-c-basic),
Copyright 2014 Melissa O'Neill, under [Apache-2.0](../LICENSES/Apache-2.0.txt).
Changes add PRISM's checked API, uniform conversion, and Box–Muller cache. That
source file retains its license notice; the other new PRISM code is MIT licensed.
No external source is fetched by the build.
