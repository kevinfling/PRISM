# PRISM

**Partitioned Reconstruction and Inference of Spectral Modes** is a C11 library
for one-dimensional, deterministic DCT filter banks. **v0.2.0-dev** adds opt-in
seeded covariance sampling, caller-provided innovation sources, an optional
pinned SIGMA normal adapter, and a fixed-linear ensemble research driver to the
v0.1 filter banks. Fixed bands and differences of heat smoothers retain
caller-selected parameters, explicit components, and a mandatory residual.
This is a development snapshot, not a published v0.2.0 release.

For those fixed parameters, the components are linear filtered views of the input.
Their weights sum to one, so they reconstruct the input in exact arithmetic.
The current tests measure how closely the implementation realizes that identity.
They do not establish physical source recovery, denoising benefit, or statistical
inference. The word *inference* in the project name describes a research direction;
the implementation estimates no latent model, probabilities, or calibrated uncertainty.

[Thesis and null hypotheses](#thesis-and-null-hypotheses) ·
[Build](#build-and-install) · [Example](#example) ·
[Mathematical specification](#mathematical-specification) ·
[Evidence](#evidence-and-release-status)

## Thesis and null hypotheses

Our engineering thesis is that a small, inspectable spectral reconstruction core
can be a useful baseline for studying multiscale signal representations. The
research question is whether future partition selection, perturbations, or
state-space models improve a **prespecified downstream task** enough to justify
their cost. That question is open. The mathematical identities below are standard
linear algebra and spectral filtering properties, with no novelty claim.

We use two distinct null hypotheses:

1. **Application null:** at matched resource budgets and on held-out data, a
   proposed method has no practically meaningful advantage over a named simpler
   baseline. If a higher-is-better score is used, define the paired improvement
   $`\Delta`$ and the minimum useful effect $`\delta\gt 0`$ before final evaluation:
   $`H_0:\mathbb E[\Delta]\le\delta`$. For a loss, reverse the sign consistently.
   The task, baseline, primary metric, compute budget, split, tuning procedure,
   and uncertainty calculation must be fixed in the experiment protocol. PRISM has
   not run such an application experiment or rejected this null.
2. **Fixed-linear perturbation null:** independent, zero-mean additive trials
   passed through the same fixed linear operator preserve its deterministic
   expectation and add estimator variance. This is a conditional mathematical
   result, derived in [§10](#10-the-fixed-linear-ensemble-identity-research-only),
   rather than an empirical discovery. The opt-in v0.2 development experiment
   checks this identity with zero, white, and fGn controls; it evaluates no
   application benefit. Nonlinear or adaptive extensions require separate tasks
   and protocols.

An inability to reject an application null would not prove equivalence. A
successful reconstruction test also cannot reject a source-separation or
application null: retaining every band plus the residual retains the input.

## Build and install

Requirements: CMake 3.16+, a C11 compiler, the C math library, and, by default,
pkg-config plus a separately installed double-precision FFTW3. On Debian/Ubuntu,
the FFTW development package is `libfftw3-dev`. Configuration downloads nothing.
Local validation covers Linux/aarch64 with GCC and Clang; the CI configuration
also targets Linux/x86-64. Remote CI results and other platforms are not yet
established by the local validation record.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel
ctest --test-dir build -L prism --output-on-failure
./build/examples/prism_example fixed
./build/examples/prism_example heat
cmake --install build --prefix "$PWD/install"
```

Use `-DPRISM_ENABLE_FFTW=OFF` to build with only the quadratic reference transform;
that configuration needs no third-party transform library or pkg-config.
The reference backend is intended for small inputs and validation.
The [SIGMA setup and audit](docs/SIGMA.md) describes its separately supplied,
hash-pinned dependency; the default build requires no SIGMA source.

| CMake option | Default | Effect |
| --- | --- | --- |
| `PRISM_ENABLE_FFTW` | ON | Build the FFTW backend and require system FFTW3 |
| `PRISM_ENABLE_SIGMA` | OFF | Enable the pinned optional normal-source adapter |
| `PRISM_SIGMA_SOURCE_DIR` | empty | Source directory with the audited SIGMA headers/license; required only when enabled |
| `PRISM_BUILD_TESTS` | ON for a standalone build | Build and register PRISM checks |
| `PRISM_BUILD_EXAMPLES` | ON for a standalone build | Build `prism_example` |
| `PRISM_BUILD_BENCHMARKS` | ON for a standalone build | Build the POSIX decomposition and innovation benchmarks |
| `PRISM_BUILD_EXPERIMENTS` | OFF | Build the POSIX fixed-linear research driver and register its checks |
| `PRISM_WARNINGS_AS_ERRORS` | OFF | Treat project warnings as errors with GCC/Clang |
| `PRISM_SANITIZE` | OFF | Instrument with ASan/UBSan using GCC/Clang |
| `BUILD_SHARED_LIBS` | OFF | Select shared rather than static library output |

Tests, examples, and benchmarks default to OFF when PRISM is included with
`add_subdirectory`. No native-CPU or fast-math flags are forced.

An installed consumer uses:

```cmake
find_package(PRISM 0.2 CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE prism::prism)
```

Point `CMAKE_PREFIX_PATH` at the installation prefix. Public headers support C++
inclusion; the implementation is C11. This is a pre-1.0 API: compatibility is
scoped to the 0.2 series, and later minor versions may change the API or ABI.

For instrumented checks:

```sh
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug \
    -DPRISM_SANITIZE=ON -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build-sanitize --parallel
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    ctest --test-dir build-sanitize -L prism --output-on-failure
```

Leak detection requires an environment permitting the sanitizer's process
inspection. See the [validation record](docs/VALIDATION.md) for local commands,
coverage, and limitations.

## Example

This complete program uses the reference backend explicitly, so it also works
when FFTW is disabled. The result includes three components **and** the remaining
low-frequency content. Capacities count doubles, not bytes.

```c
#include <math.h>
#include <prism.h>
#include <stdio.h>

int main(void) {
    enum { N = 257, J = 3 };
    size_t edges[] = {N, 96, 24, 4};
    prism_config config = {
        .n = N,
        .component_count = J,
        .sample_spacing = 1.0,
        .mode = PRISM_BANDS_FIXED,
        .band_edges = edges
    };
    double x[N], components[J * N], residual[N], reconstructed[N];
    for (size_t i = 0; i < N; ++i)
        x[i] = 2.0 + sin(0.4 * (double)i) + 0.3 * cos(1.7 * (double)i);
    prism_result result = {
        .components = components,
        .component_capacity = J * N,
        .residual = residual,
        .residual_capacity = N
    };
    prism_plan *plan = NULL;
    prism_workspace *workspace = NULL;
    int status = prism_plan_create(&config, prism_backend_reference(), NULL, &plan);
    if (status == PRISM_OK)
        status = prism_workspace_create(plan, &workspace);
    if (status == PRISM_OK)
        status = prism_decompose(plan, workspace, x, &result);
    if (status == PRISM_OK)
        status = prism_reconstruct(&result, reconstructed, N);
    if (status == PRISM_OK)
        printf("relative reconstruction error: %.3e\n",
               result.reconstruction_relative_error);
    else
        fprintf(stderr, "PRISM: %s\n", prism_status_string(status));
    prism_workspace_destroy(workspace);
    prism_plan_destroy(plan);
    return status == PRISM_OK ? 0 : 1;
}
```

To select heat partitions, set `.mode = PRISM_BANDS_HEAT_TELESCOPE` and supply
`.heat_times` pointing to J increasing positive diffusion times, for example
`{0.5, 8.0, 128.0}`. The unused band-edge pointer may be NULL. Setting the backend
to NULL selects FFTW when built in, otherwise reference.

Plans own their precomputed weights; configuration arrays can be released after
creation. The plan must outlive its workspaces. Each simultaneous execution needs
its own workspace and output buffers. PRISM's execution paths perform no explicit
allocation; custom backends must satisfy the same contract. FFTW planning and
destruction require caller serialization. Lifetime, buffer, callback, and error
contracts are specified in [docs/API.md](docs/API.md).

## Mathematical specification

The following treatment specifies the complete implemented finite-dimensional
operators, their reconstruction and stability properties, and their numerical
limits. Unless explicitly discussing floating point, equations are over the real
numbers. Indices for samples and spectral coefficients start at zero; there are
J components with $`j=0,\ldots,J-1`$, followed by one residual.

### 1. Samples, units, and the orthonormal DCT

Let $`n\ge2`$, $`d_x\gt 0`$, $`L=nd_x`$, and $`x\in\mathbb R^n`$. Sample locations are
$`s_i=(i+\tfrac12)d_x`$ for $`i=0,\ldots,n-1`$. Define

```math
C_{ki}=a_k\cos\!\left(\frac{\pi k(i+\tfrac12)}{n}\right),\qquad
 a_0=\frac1{\sqrt n},\quad a_k=\sqrt{\frac2n}\ (k\gt 0).
```

Analysis and synthesis are

```math
\widehat x=Cx,\qquad x=C^{\mathsf T}\widehat x.
```

To verify normalization, the half-grid cosine sum satisfies

```math
\sum_{i=0}^{n-1}\cos\!\left(\frac{\pi m(i+\tfrac12)}n\right)=0,
\qquad m=1,\ldots,2n-1.
```

This follows from a finite geometric series, or the cosine-sum identity
$`\sin(n\theta)/(2\sin(\theta/2))`$ with $`\theta=\pi m/n`$.
The product-to-sum formula then gives unscaled cosine inner products n at
$`k=\ell=0`$, n/2 at $`k=\ell\gt 0`$, and zero otherwise. Hence
$`CC^{\mathsf T}=C^{\mathsf T}C=I`$ and

```math
\|x\|_2^2=\sum_{k=0}^{n-1}|\widehat x_k|^2.
```

DCT types and boundary placements are classical; see
[Strang, *The Discrete Cosine Transform*](https://math.mit.edu/~gs/papers/dct.pdf).

### 2. Boundary convention and continuum interpretation

The unique cosine polynomial in the retained n-dimensional space interpolating
the data is

```math
f_x(s)=\sum_{k=0}^{n-1} a_k\widehat x_k\cos(\pi ks/L).
```

It satisfies $`f_x(s_i)=x_i`$ and $`f_x'(0)=f_x'(L)=0`$. This is a modeling choice
about the interpolant and even reflection at the interval boundaries; arbitrary
input samples are accepted. It does not show that the unknown process generating
the samples has zero boundary derivatives or that the data resolve its spectrum.
The physical wavenumber is $`\omega_k=\pi k/L`$, and its frequency in cycles per
unit length is $`k/(2L)`$. For temporal signals, length units here are time units.

Writing continuous orthonormal cosine functions as
$`\phi_0=L^{-1/2}`$ and $`\phi_k(s)=\sqrt{2/L}\cos(\pi ks/L)`$ gives coefficients
$`\sqrt{d_x}\,\widehat x_k`$ and

```math
\int_0^L |f_x(s)|^2\,ds=d_x\|x\|_2^2.
```

Thus Euclidean energy and physical integral energy differ by $`d_x`$ under this
interpolation convention. Global cosine filters use the whole record, generally
have boundary effects, and are not causal streaming estimators or exactly
translation-equivariant filters on the finite interval.

### 3. General partitions and reconstruction

Let $`w_0,\ldots,w_{J-1},w_R\in\mathbb R^n`$ obey, for every k,

```math
w_R[k]+\sum_{j=0}^{J-1}w_j[k]=1.
```

For each label $`a\in\{0,\ldots,J-1,R\}`$, set

```math
A_a=C^{\mathsf T}\mathop{\mathrm{diag}}(w_a)C,\qquad
c_j=A_jx,\qquad r=A_Rx.
```

Summing gives

```math
\sum_j c_j+r
=C^{\mathsf T}\mathop{\mathrm{diag}}\!\left(\sum_jw_j+w_R\right)Cx
=C^{\mathsf T}Cx=x.
```

This identity holds for every x and every valid partition in exact arithmetic.
It establishes reconstruction, with no source-identification conclusion.
Here *component* means this particular filtered output; *residual* means the
retained complement. The residual need not be noise or an error term.

### 4. Fixed nonoverlapping bands

The configuration supplies integer edges

```math
n=b_0\gt b_1\gt \cdots\gt b_J\ge1.
```

Weights are

```math
w_j[k]=\mathbf1_{\{b_{j+1}\le k\lt b_j\}},\qquad
w_R[k]=\mathbf1_{\{0\le k\lt b_J\}}.
```

These disjoint sets cover all n indices, proving the partition property.
Components are ordered from high to low frequency and DC remains in the residual.
Because weights are binary and disjoint,

```math
A_a^{\mathsf T}=A_a,\quad A_a^2=A_a,\quad A_aA_b=0\quad(a\ne b).
```

They are orthogonal projections onto the selected DCT subspaces. In particular,
$`\|x\|_2^2=\sum_j\|c_j\|_2^2+\|r\|_2^2`$. Abrupt spectral cutoffs can produce
ringing near discontinuities. Edge indices are user parameters, not inferred
physical mode boundaries; the API requires $`1\le J\le n-1`$ implicitly through
strictly descending integer edges ending at least at 1.

### 5. Telescoping heat partitions

For $`0=t_0\lt t_1\lt \cdots\lt t_J`$, define the finite spectral operator

```math
A=C^{\mathsf T}\mathop{\mathrm{diag}}(\lambda_k)C,\quad
\lambda_k=(\pi k/L)^2,\quad
S(t)=e^{-tA}=C^{\mathsf T}\mathop{\mathrm{diag}}(e^{-t\lambda_k})C.
```

On the cosine polynomial space, this is the Neumann heat evolution
$`\partial_t f=\partial_s^2 f`$ with initial interpolant $`f_x`$. Its finite spectral
representation is exact in real arithmetic; no accuracy claim follows for an
unresolved underlying continuous signal. Diffusion time t has units of length
squared with diffusivity set to one. For diffusivity $`\kappa`$, use $`t=\kappa\tau`$.

The implemented weights and outputs are

```math
w_j[k]=e^{-t_j\lambda_k}-e^{-t_{j+1}\lambda_k},\quad
w_R[k]=e^{-t_J\lambda_k},
```

```math
c_j=(S(t_j)-S(t_{j+1}))x,\qquad r=S(t_J)x.
```

Their sum telescopes from $`S(0)=I`$, proving reconstruction. All weights lie in
$`[0,1]`$, component DC weights are zero, and $`w_R[0]=1`$. Semigroup and limiting
properties follow coefficient by coefficient:

```math
S(t)S(u)=S(t+u),\quad S(0)=I,\quad
\lim_{t\to\infty}S(t)x=\overline x\,\boldsymbol1.
```

The first band is a high-pass complement. For a later band with endpoints
$`0\lt a\lt b`$, its continuous-in-$`\lambda`$ response has its maximum at

```math
\lambda_* = \frac{\log(b/a)}{b-a}.
```

This follows by differentiating $`e^{-a\lambda}-e^{-b\lambda}`$; the actual response
is sampled only at the retained $`\lambda_k`$. Adjacent responses overlap and are
generally neither projections nor mutually orthogonal. Nonnegative spectral
weights imply positive semidefinite matrices, not nonnegative matrix entries,
pointwise-positive outputs, or a discrete maximum principle. Smooth spectral
attenuation alone does not establish recovery of a slowly varying source.

For small diffusion increments, code evaluates the difference as

```math
e^{-t_j\lambda_k}\left[-\mathop{\mathrm{expm1}}(-(t_{j+1}-t_j)\lambda_k)\right]
```

rather than subtracting nearly equal exponentials. If $`d_x`$ is multiplied by q
and all times by $`q^2`$, the ideal weights are unchanged. Fixed-band weights are
independent of spacing; spacing only gives their physical frequency interpretation.

### 6. Continuum spectral versus finite-difference heat

The cell-centered Neumann second-difference matrix has diagonal
$`(1,2,\ldots,2,1)/d_x^2`$ and off-diagonals $`-1/d_x^2`$. The same C diagonalizes it,
but its eigenvalues are

```math
\lambda_k^{\mathrm{FD}}=\frac4{d_x^2}\sin^2\!\left(\frac{\pi k}{2n}\right).
```

For $`k\gt 0`$,

```math
\frac{\lambda_k^{\mathrm{FD}}}{\lambda_k}
=\left(\frac{\sin z_k}{z_k}\right)^2,\qquad z_k=\frac{\pi k}{2n}.
```

The ratio tends to one at low relative frequencies, but generally differs near
the top of the spectrum. PRISM uses $`\lambda_k`$, **not**
$`\lambda_k^{\mathrm{FD}}`$. Sharing a cosine basis does not make these two heat
operators identical. Neither boundary model is automatically suitable for every
data-generating process.

### 7. Energy, stability, and what a component can mean

Let $`y_a=A_ax`$, including the residual. For real weights, orthonormality gives

```math
\langle y_a,y_b\rangle=\sum_k w_a[k]w_b[k]|\widehat x_k|^2,
```

```math
\|x\|_2^2=\sum_a\|y_a\|_2^2+2\sum_{a\lt b}\langle y_a,y_b\rangle.
```

For the implemented nonnegative partitions, these cross terms are nonnegative.
They vanish for disjoint fixed bands; heat-band energies generally sum to less
than input energy. Reporting them as independent fractions of explained energy
would therefore be misleading. Orthogonality also does not imply statistical
independence without additional distributional assumptions.

For either implemented family with its parameters held fixed, the operator norm is

```math
\|A_a\|_{2\to2}=\max_k|w_a[k]|\le1,
\qquad \|A_a(x+e)-A_ax\|_2\le\|e\|_2.
```

This is absolute perturbation stability of each linear filter in exact
arithmetic. Relative error can still be large for a component of nearly zero
energy. For the stacked representation $`Fx=(y_0,\ldots,y_{J-1},r)`$,

```math
\|Fx\|_2^2=\sum_k q_k|\widehat x_k|^2,\qquad q_k=\sum_a w_a[k]^2,
\qquad \frac1{J+1}\le q_k\le1.
```

The lower bound is Cauchy–Schwarz applied to J+1 weights summing to one; the upper
bound uses nonnegativity. Thus the representation is injective and has bounded
norms, but is not generally a Parseval or tight frame. These elementary
properties are distinct from a measured stability benefit on real data.

Linearity and DC placement give, for fixed parameters,

```math
c_j(\alpha x+\beta\boldsymbol1)=\alpha c_j(x),\qquad
r(\alpha x+\beta\boldsymbol1)=\alpha r(x)+\beta\boldsymbol1.
```

Parameter selection based on x would generally make the overall procedure
nonlinear; PRISM performs no such selection. If $`x=s_1+s_2`$, then
$`x=(s_1+h)+(s_2-h)`$ for any h, so observations alone cannot identify the two
sources. A spectral partition supplies a representation choice, not the missing
source assumptions. Global bands cannot in general separate sources occupying
the same spectral region, including crossing or overlapping chirps.

### 8. Floating-point reconstruction and diagnostic limits

For stored weights, define partition defect

```math
\delta_k=\sum_a \widetilde w_a[k]-1,\qquad d=\max_k|\delta_k|.
```

With otherwise exact transforms, the reconstructed-minus-input vector is
$`C^{\mathsf T}\mathop{\mathrm{diag}}(\delta)Cx`$, so

```math
\|\widetilde x-x\|_2\le d\|x\|_2.
```

To separate numerical contributions, write computed analysis as $`Cx+e_f`$ and
each computed synthesis as
$`C^{\mathsf T}\mathop{\mathrm{diag}}(\widetilde w_a)(Cx+e_f)+e_a`$.
Including final accumulation error $`e_\Sigma`$ gives

```math
\|\widetilde x-x\|_2
\le d\|x\|_2+(1+d)\|e_f\|_2+\sum_a\|e_a\|_2+\|e_\Sigma\|_2.
```

Here $`e_a`$ absorbs multiplication and synthesis errors. This accounting identity
and bound do not provide a machine-independent bound on the actual error terms.
Underflow, overflow, transform normalization, and cancellation matter.

The implementation reports

```math
\rho=\begin{cases}\|x-\widetilde x\|_2/\|x\|_2,&\|x\|_2\gt 0,\\ \|\widetilde x\|_2,&x=0,\end{cases}
\qquad m=\max_i|x_i-\widetilde x_i|.
```

It uses long-double accumulation and hypot norms; long double need not have
extra precision on every platform. Diagnostics are numerical estimates, not
certified error bounds. They use an extended-precision sum of the stored outputs;
`prism_reconstruct` additionally rounds that sum into its double output buffer.

The test targets are $`\rho\le10^{-10}`$ on the documented test domain and partition
defect below $`10^{-14}`$. **PRISM_OK does not enforce these tolerances.** It means
the checked execution completed with finite outputs and diagnostics. Callers may
inspect the diagnostics and apply their own acceptance criteria. Extreme finite
inputs can fail even when an ideal normalized operator is bounded: intermediate
unnormalized transforms and norm calculations also need representable values.
No bitwise reproducibility across compilers, libm implementations, planner state,
backends, or platforms is promised.

### 9. Backend normalization and computational cost

The reference backend evaluates cosine sums using periodically re-anchored
trigonometric rotation and compensated accumulation. The test oracle separately
evaluates a direct cosine matrix. This supplies an implementation cross-check;
it is not third-party validation or a formal verification of either program.

FFTW's DCT-II (`REDFT10`) yields
$`Y_k=2\sum_i x_i\cos(\pi k(i+\tfrac12)/n)`$.
PRISM obtains $`\widehat x_0=Y_0/(2\sqrt n)`$ and
$`\widehat x_k=Y_k/\sqrt{2n}`$ for $`k\gt 0`$.
For DCT-III (`REDFT01`), it supplies
$`Z_0=\widehat x_0/\sqrt n`$ and $`Z_k=\widehat x_k/\sqrt{2n}`$, producing
$`Z_0+2\sum_{k\gt 0}Z_k\cos(\pi k(i+\tfrac12)/n)=x_i`$.
These factors adapt FFTW's unnormalized pair to C and its transpose; see the
[FFTW transform definitions](https://www.fftw.org/fftw3_doc/1d-Real_002deven-DFTs-_0028DCTs_0029.html).

One decomposition performs one analysis and J+1 syntheses, plus $`O(Jn)`$ work for
weighting, validation, and diagnostics. If one transform costs $`D(n)`$, execution
cost is $`(J+2)D(n)+O(Jn)`$. The reference implementation has $`D(n)=O(n^2)`$;
FFTW's real-to-real implementation uses $`O(n\log n)`$ algorithms, including for
prime sizes ([FFTW manual](https://www.fftw.org/fftw3_doc/More-DFTs-of-Real-Data.html)).
These are operation-count descriptions, not measured latency guarantees.
Weight preparation costs $`O(Jn)`$ plus backend planning at workspace creation.

A plan stores (J+1)n double weights. A workspace stores two n-double arrays plus
backend state; the FFTW backend adds one n-double buffer and its plans. Caller
components and residual use another (J+1)n doubles. Output construction is
allocation-free within PRISM after planning; a caller-supplied backend must also
honor that contract. The FFTW adapter uses `ESTIMATE | UNALIGNED`, preserving
support for ordinarily aligned allocator memory. Planner state can affect the
selected algorithm. See [FFTW planner flags](https://www.fftw.org/fftw3_doc/Planner-Flags.html)
and [thread-safety requirements](https://www.fftw.org/fftw3_doc/Thread-safety.html).

### 10. The fixed-linear ensemble identity (research only)

This identity is checked by the optional [research driver](experiments/ensemble.c);
it is not an ensemble mode of `prism_decompose`.
Let S be a fixed linear smoother, q a fixed nonnegative integer,
$`T=(I-S)^q`$ (with $`T=I`$ when q=0), r a fixed input, and $`\eta_e`$ independent random vectors with mean
zero and common finite covariance $`\Sigma`$. For a noise amplitude $`\sigma`$
independent of ensemble count E, define

```math
z_E=\frac1E\sum_{e=1}^E T(r+\sigma\eta_e)
=Tr+\sigma T\overline\eta_E.
```

By linearity and independence,

```math
\mathbb E[z_E]=Tr,\qquad
\mathop{\mathrm{Cov}}(z_E)=\frac{\sigma^2}E T\Sigma T^{\mathsf T}.
```

Gaussianity is unnecessary for these first two moments. For a fixed target u,
expanding squared error and using the zero cross expectation gives

```math
\mathbb E\|z_E-u\|_2^2
=\|Tr-u\|_2^2+\frac{\sigma^2}E\mathop{\mathrm{tr}}(T\Sigma T^{\mathsf T}).
```

Thus this particular independent perturbation procedure cannot reduce expected
squared error relative to Tr for that fixed r, u, and T. If each trial instead
uses amplitude $`\sigma/\sqrt E`$, the covariance term scales as $`E^{-2}`$.
Correlated or antithetic trials require their cross-covariances; zero-mean noise
alone does not establish the independent-trial formula.

If S is DCT-diagonal with multiplier $`s_k`$, T has multiplier $`(1-s_k)^q`$ and can
be applied deterministically. Adaptive stopping, input-dependent filters,
nonlinear estimators, or a different downstream objective fall outside this
comparison. This identity does not rule out benefits for those procedures; they
need separately specified experiments. Nor does DCT-filtering white noise
produce an arbitrary target Toeplitz covariance: its covariance has the form
$`C^{\mathsf T}\mathop{\mathrm{diag}}(v_k)C`$ for chosen spectral variances $`v_k`$.
The opt-in sampler below constructs the specified Toeplitz target through a
checked covariance factorization or circulant embedding.

### 11. Target covariance and Gaussian sampling

The opt-in [noise API](include/prism/noise.h) accepts stationary target lags
$`\gamma(0),\ldots,\gamma(n-1)`$ and defines the Toeplitz matrix
$`K_{ij}=\gamma(|i-j|)`$. For unit-spaced fractional Gaussian increments,

```math
\gamma_H(0)=v,\qquad
\gamma_H(k)=\frac v2\left[(k+1)^{2H}-2k^{2H}+(k-1)^{2H}\right],\quad k\ge1,
\qquad 0\lt H\lt1,\ v\ge0.
```

Here v is increment variance, not a spectrum normalization inferred from H.
At H=1/2, K=vI. An fBm sampling interval h with amplitude a requires explicitly
choosing $`v=a^2h^{2H}`$. The implementation returns increments, not fBm levels.
The [sampling specification](docs/STOCHASTIC.md#covariance-definition-and-scale)
derives the stable even-binomial evaluation used instead of cancellation-prone
power subtraction.

For dense Cholesky, construct $`K=LL^{\mathsf T}`$ and return $`Lz`$ for ideal
$`z\sim N(0,I_n)`$. Computed negative pivots are rejected. Zero pivots require
zero corresponding remainders; no implicit jitter is added. Setup costs O(n³),
sampling O(n²), and the stored factor n² doubles.

For circulant embedding, set m=2(n−1), or m=1 when n=1, and mirror the lags:
$`c=(\gamma(0),\ldots,\gamma(n-1),\gamma(n-2),\ldots,\gamma(1))`$.
For n=1 use only the lag-zero entry; n=2 has no appended mirror tail.
With the unnormalized complex DFT matrix
$`F_{kj}=e^{-2\pi\mathrm{i}kj/m}`$, eigenvalues $`\lambda=Fc`$, and
$`\lambda_k\ge0`$, define

```math
G=\frac1m F^*\mathop{\mathrm{diag}}(\lambda)F,\qquad
B=\frac1m F^*\mathop{\mathrm{diag}}(\sqrt\lambda)F.
```

Since $`FF^*=mI`$, $`BB^{\mathsf T}=G`$. For ideal independent Gaussian z of
length m, return the first n entries of Bz. The corresponding block of G is K,
so the returned covariance is K in exact arithmetic. Forward DFT, square-root
weighting, and inverse DFT divided by m implement B. This uses a circulant
embedding of the target, rather than treating the DCT as a Toeplitz diagonalizer.
The reference DFT costs O(m²); optional FFTW supplies O(m log m) transform
algorithms. Plans and workspaces use O(m) PRISM storage for this path.

A valid Toeplitz target need not have a PSD **chosen embedding**. For example,
$`\gamma=(1.05,0.9,0.62)`$ has a positive definite 3×3 target, while its mirror
embedding has eigenvalues 3.47, 0.43, −0.13, 0.43. Embedding rejection therefore
does not establish invalidity of the target. There is no automatic enlargement
or fallback.

The default eigenvalue tolerance is zero. An explicit caller tolerance τ in
[0,1) permits a negative computed eigenvalue to be rounded to zero only when
$`-\lambda_k\le\tau\max_j|\lambda_j|`$. Diagnostics report the minimum/maximum
computed eigenvalues, imaginary roundoff, corrected-bin count with conjugate
multiplicities, and maximum correction. Such rounding changes the target:

```math
\Delta G=\frac1mF^*\mathop{\mathrm{diag}}(\delta\lambda)F,\qquad
\|\Delta G\|_2=\max_k|\delta\lambda_k|.
```

Cropping cannot increase this exact-arithmetic spectral-norm bound. It does not
include other floating-point errors. The copied target lags remain the requested
lags, not the corrected covariance. A corrected plan is not an exact-target
sampler. See [full factorization and error contracts](docs/STOCHASTIC.md).

### 12. Seeded execution and statistical checks

A caller-owned [PCG RNG](include/prism/rng.h) takes an explicit seed and logical
stream ID. Two sequential 32-bit words produce a 53-bit [0,1) uniform. For two
such uniforms U,V, Box–Muller returns

```math
R=\sqrt{-2\log(1-U)},\qquad \Theta=2\pi V,\qquad
(Z_1,Z_2)=(R\cos\Theta,R\sin\Theta).
```

The second value is cached; chunking draws preserves the sequence. Integer
outputs are specified exactly; normals depend on libm. Separate RNGs and
workspaces permit concurrent sampling with one immutable plan. Logical stream
IDs make scheduling reproducible; different IDs alone do not prove independence.
There is no clock seeding or global RNG. Finite-grid pseudorandom normals and
floating-point factorizations approximate the ideal Gaussian construction;
“exact covariance” does not promise literal continuous Gaussian samples or
bitwise portability. Sampling allocates no PRISM memory after setup.

For M ideal independent, zero-mean Gaussian draws, define the known-zero-mean
second-moment estimate $`\widehat K=M^{-1}\sum_{s=1}^M y_sy_s^{\mathsf T}`$.
Gaussian fourth moments give

```math
\mathbb E\widehat K=K,\qquad
\mathop{\mathrm{Var}}(\widehat K_{ij})=
\frac{K_{ii}K_{jj}+K_{ij}^2}{M},\qquad
\mathbb E\|\widehat K-K\|_F^2=
\frac{(\mathop{\mathrm{tr}}K)^2+\mathop{\mathrm{tr}}(K^2)}{M}.
```

Tests use those scales for selected-entry six-standard-error checks and a
separate broad full-matrix second-moment check. They also compare an independent
full cosine-matrix circulant square root, published PCG integer vectors,
allocation failures, capacities, and serial/concurrent logical-stream execution.
Finite fixed-seed checks do not prove the assumed independence or exactness.

The optional experiment uses q=3, heat time 16, dx=1, n=64, fixed σ=0.25,
20 seeded replicates, E=1,4,16,64,256, and zero/white/fGn controls. With
$`A=TKT^{\mathsf T}`$, ideal Gaussian noise gives

```math
D_E=\|z_E-Tx\|_2^2,\qquad
\mathbb E D_E=\frac{\sigma^2}E\mathop{\mathrm{tr}}A,\qquad
\mathop{\mathrm{Var}}(D_E)=\frac{2\sigma^4}{E^2}\mathop{\mathrm{tr}}(A^2).
```

This follows by diagonalizing A and using variance two for squared independent
standard normals. A direct cosine matrix supplies the prediction independently
of the transform implementation. The [fixed protocol](docs/STOCHASTIC.md#fixed-linear-protocol-and-interpretation)
specifies input, streams, nested-prefix coupling, raw CSV fields, timings, and
uncertainty limits. Enable it with `-DPRISM_BUILD_EXPERIMENTS=ON` and run
`build/experiments/prism_ensemble_experiment fftw` or `reference`.

The report's approximate t intervals describe variation across seeded
replicates, assuming effective independence. They do not calibrate physical
source uncertainty. Reduced perturbation error with larger E compares noisy
estimates with each other; deterministic Tx already has zero perturbation error.
This experiment supplies no application evaluation or ensemble advantage over Tx.

### 13. Innovation boundary and optional SIGMA source

The [innovation interface](include/prism/innovation.h) separates a caller-owned
normal source from the covariance transformation. A descriptor supplies state,
its byte extent, and a fill callback; `prism_rng_source` adapts PCG.
`prism_noise_sample_from` stages a source's innovations before producing a
sample. `prism_noise_transform` instead applies the square root to explicit
caller-provided innovations, with no RNG consumption. Required innovation count
is n for dense Cholesky and the embedding length m for a circulant plan.

For any supplied vector z with mean μ and covariance V and the square-root
operator R (including cropping when needed),

```math
mathbb E[Rz]=Rmu,qquad mathop{mathrm{Cov}}(Rz)=RVR^{mathsf T}.
```

For an uncorrected exact target factor, centered, identity-covariance innovations
give the target second moments; Gaussian innovations additionally give the ideal
Gaussian field. Optional eigenvalue correction changes that target as in §11. Finite-value
validation cannot verify those distributional assumptions for an arbitrary
callback. Source, state, and output must be disjoint, callbacks must not allocate,
and no sampling-time planning occurs. Callback failures may consume source
state; covariance sampling still leaves the caller's sample unchanged.

The opt-in [SIGMA adapter](include/prism/sigma.h) pins revision
`6e9a1ef526365474738d80d4eb71fd945a13abd6` and checks consumed-header/license
hashes. It uses caller-owned xoshiro256** state and scalar cached polar normals:
for independent ideal uniforms, form u=2U−1, v=2V−1, reject until 0<a=u²+v²<1,
and return

```math
(Z_1,Z_2)=(u,v)sqrt{-2log(a)/a}.
```

One value is cached. The adapter fills buffers by scalar calls so chunking
preserves that sequence, including cached draws. It exposes no SIGMA types or
headers to consumers and imports no dependency optimization flags.
The [complete seed/state specification and upstream audit](docs/SIGMA.md)
explain why upstream buffered/SIMD sequence equivalence is not assumed.
Different stream keys do not prove independent random variables. Assign keys
to logical samples for scheduling invariance, and separately define aggregation
order if requiring bitwise repeatability. No cross-platform normal-value
bitwise identity is promised.

Tests compare word transitions with an independent oracle, check scalar/chunked
normals, raw moments and fixed-grid CDF values, and compare serial scheduling
with 1,2,4,7 threads. At lengths 257,511,513, all innovation basis vectors recover
the FFTW sampler's numerical square root. Boundary-row and diagonal inner
products are checked against dense Toeplitz entries before Monte Carlo
lag-covariance checks. Those 4,096-draw checks report approximate pointwise 95%
intervals from the Gaussian fourth-moment variance in §12; they are not
simultaneous confidence bands. The reference-only build additionally checks a
dense n=513, H=0.7 factor column. These finite checks supply no all-size,
independence, or calibrated application-uncertainty guarantee.

The [innovation benchmark](benchmarks/bench_innovations.c) separates normal-fill
cost from correlated-noise sampling plus perturbation and full decomposition.
Measured throughput is not evidence of better decomposition. The optional
adapter uses scalar normal draws; its measurements are not SIMD normal results.
Other distributions, online-statistics adapters, antithetic/common-random-number
comparisons, and auxiliary operators remain unimplemented.

## Evidence and release status

| Statement | Basis | Boundary of the claim |
| --- | --- | --- |
| Components plus residual reconstruct x | Exact-arithmetic derivation | Floating-point behavior is measured separately |
| Each fixed filter is nonexpansive in L2 | Operator-norm calculation | Does not bound relative error of a tiny component or prove source separation |
| Tests meet the chosen reconstruction tolerance | Local test observations | Finite test domain, compilers, and machine; no all-input guarantee |
| Capacity/overlap checks and sanitizer checks pass | Source checks and executed tests | Not a formal safety proof or exhaustive race analysis |
| Same parameters specify the same real operator | Mathematical definition | No cross-platform bitwise reproducibility claim |
| Denoising or application advantage | Open hypothesis | No held-out application evaluation |
| Covariance and fixed-linear ensemble checks | Specified construction, internal oracles, seeded simulation | Floating-point and PRNG limitations; no application or all-input guarantee |
| Comparative speed or new mathematics | Not established or claimed | Smoke timings and elementary identities are insufficient evidence |

This development snapshot implements the specified filter banks and opt-in
covariance samplers. It has no adaptive partition selection, local frames,
fractional calculus, cepstrum, RTS model, EMD/envelope sifting, production ensemble
decomposition, or calibrated uncertainty quantification API. Reconstruction and component energy do not demonstrate
separation quality, denoising benefit, or physical interpretability.

The [validation record](docs/VALIDATION.md) distinguishes test coverage and local
observations from untested platforms and research hypotheses. Tests cover both
partitions and backends at lengths 128,257,511,512,513,1024,4096,16384, alongside
small-length oracles, seeded properties, capacity/allocation failures, and
concurrent workspaces. Passing those tests does not establish all-input
correctness, universal numerical conditioning, or production readiness.

The [benchmark](benchmarks/bench_prism.c) separates setup, first execution, and
warm execution including diagnostics. Use the same size and repetition count for
backend comparisons, for example:

```sh
./build/benchmarks/bench_prism fftw 512 20
./build/benchmarks/bench_prism reference 512 20
```

The archived [smoke results](results/2026-10-21-jetson-orin-nano/README.md) use
unequal sizes and one uncontrolled session. They support no relative-speed,
optimality, or cross-platform performance claim. New performance claims require
matched configurations, repeated sessions, raw measurements, and uncertainty.

[CONTRIBUTING.md](CONTRIBUTING.md) describes checks and evidence expectations;
[CHANGELOG.md](CHANGELOG.md) distinguishes initial release history from current
development. Other backend adapters, auxiliary operators, and adaptive models
remain later work.

## License and attribution

PRISM's code is [MIT licensed](LICENSE), Copyright (c) 2026 Kevin Fling, except
`src/rng.c`, which adapts Melissa O'Neill's PCG code under
[Apache-2.0](LICENSES/Apache-2.0.txt), Copyright 2014 Melissa O'Neill. See
[attribution](docs/STOCHASTIC.md#attribution).
DCT filter banks, heat semigroups, and the reconstruction identity are established
techniques; the references above provide context for the normalization and
boundary conventions. No original mathematical result is claimed here.

SIGMA is a separately supplied optional MIT dependency, Copyright (c) 2025
Kevin Fling. Its license is installed with enabled builds; the
[adapter audit](docs/SIGMA.md) records the pin and supported scalar API.

FFTW is an optional, separately installed dependency. Its default distribution
uses **GPL version 2 or later**, and alternative licenses are available from its
licensor ([FFTW license](https://www.fftw.org/fftw3_doc/License-and-Copyright.html)).
PRISM's MIT license does not replace FFTW's terms for FFTW-linked distributions.
The reference-only build avoids that dependency. FFTW is by Matteo Frigo and
Steven G. Johnson; see their [FFTW3 paper](https://www.fftw.org/fftw-paper-ieee.pdf).
