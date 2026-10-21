# Optional SIGMA normal-source adapter

PRISM v0.2.0 development can use a separately supplied, pinned SIGMA RNG.
The default remains the built-in PCG source. SIGMA is not needed to build,
install, or consume the default library. This adapter provides scalar xoshiro
words and standard-normal innovations; it supplies no new decomposition method,
probability-distribution API, SIMD normal-fill acceleration, or application result.

## Dependency boundary and pin

The audited revision is `6e9a1ef526365474738d80d4eb71fd945a13abd6` of
[StraylightRunMOO/SIGMA](https://github.com/StraylightRunMOO/SIGMA).
[PrismSigma.cmake](../cmake/PrismSigma.cmake) locks the two headers used by the
adapter and the license by SHA-256. A Git checkout is convenient but not required:
a supplied source directory with matching files also works. No dependency is
fetched by configuration or compilation. A changed consumed header is rejected,
and CMake watches those files for reconfiguration.

```sh
git clone https://github.com/StraylightRunMOO/SIGMA.git ../SIGMA
git -C ../SIGMA checkout 6e9a1ef526365474738d80d4eb71fd945a13abd6
cmake -S . -B build-sigma -DCMAKE_BUILD_TYPE=Release \
    -DPRISM_ENABLE_SIGMA=ON -DPRISM_SIGMA_SOURCE_DIR="$PWD/../SIGMA" \
    -DPRISM_BUILD_EXPERIMENTS=ON -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build-sigma --parallel 3
ctest --test-dir build-sigma -L prism --output-on-failure
./build-sigma/experiments/prism_ensemble_experiment fftw sigma > ensemble-sigma.csv 2> ensemble-sigma.txt
python3 tools/summarize_ensemble.py ensemble-sigma.csv
./build-sigma/benchmarks/bench_innovations pcg 257 100
./build-sigma/benchmarks/bench_innovations sigma 257 100
```

Use `reference` for the experiment if FFTW is disabled. SIGMA's CMake project is
not imported: it sets directory-wide fast-math and ISA options that PRISM does
not require. Only private system include paths for the two consumed headers are
added. No SIGMA type, header, or link dependency leaks into the installed public
API. Its MIT license, Copyright (c) 2025 Kevin Fling, is installed as
`share/licenses/PRISM/SIGMA-MIT.txt` when the adapter is enabled.

[prism/sigma.h](../include/prism/sigma.h) is available in every build.
`prism_sigma_available()` and `prism_sigma_revision()` report availability.
Disabled seed/draw operations return `PRISM_EBACKEND`, without drawing values.
The source descriptor borrows a caller-owned `prism_sigma_rng`; separate
concurrent calls require distinct states. No source/workspace allocation occurs
while filling normals. The [innovation interface](../include/prism/innovation.h)
allows the same covariance plan to use this or another source.

## Seed, stream, and sample policy

All integer arithmetic is unsigned modulo 2⁶⁴. Define SplitMix64's output
permutation from an input x by

```math
z_0=x+\mathtt{9e3779b97f4a7c15},\quad
z_1=(z_0\mathbin{\mathrm{xor}}(z_0\gg30))\mathtt{bf58476d1ce4e5b9},
```

```math
z_2=(z_1\mathbin{\mathrm{xor}}(z_1\gg27))\mathtt{94d049bb133111eb},\qquad
M(x)=z_2\mathbin{\mathrm{xor}}(z_2\gg31).
```

Constants above are hexadecimal; $`\gg`$ denotes logical right shift. The adapter
sets the effective seed to $`s_{\mathrm{eff}}=s\mathbin{\mathrm{xor}}M(j)`$ for
64-bit seed s and stream ID j. SIGMA then fills four xoshiro words from four
successive SplitMix64 outputs starting from that effective seed. For a fixed s,
distinct j give distinct effective seeds because M is a permutation. Different
(seed,stream) pairs can collide across seeds: this mapping is not an injective
encoding of 128 bits into 64 bits. It promises neither jump-ahead subsequences
nor probabilistic independence.

The adapter uses SIGMA's xoshiro256** transition and output:

```math
r=9\,\mathop{\mathrm{rotl}}_{64}(5s_1,7),\qquad
s'_0=s_0\mathbin{\mathrm{xor}}s_1\mathbin{\mathrm{xor}}s_3,
```

```math
s'_1=s_1\mathbin{\mathrm{xor}}s_2\mathbin{\mathrm{xor}}s_0,\quad
s'_2=s_2\mathbin{\mathrm{xor}}s_0\mathbin{\mathrm{xor}}(s_1\ll17),\quad
s'_3=\mathop{\mathrm{rotl}}_{64}(s_3\mathbin{\mathrm{xor}}s_1,45).
```

Here rotl is a 64-bit rotate, and the output uses the old state. See the
[generator authors' reference](https://prng.di.unimi.it/xoshiro256starstar.c).
The integer sequence is checked against a separately written transition and
seeding oracle. Floating-point normals still depend on libm and compiler choices.

For sample-order invariance, initialize each logical sample with the same
(seed,stream) key regardless of its worker assignment. The concurrency test uses
stream ID equal to sample ID and compares serial execution with 1,2,4,7 worker
threads for both PCG and SIGMA and every enabled covariance method. Applications
with multiple experiment dimensions must define their own stable key allocation;
the adapter does not allocate unique IDs or avoid key reuse on their behalf.
Within one continuing source, normal draw order remains part of the sequence.
Parallel floating-point aggregation also requires a defined reduction order if
bitwise equality is required; source indexing alone does not enforce it.

## Scalar and buffered policy from the audit

At the pin, the relevant signatures are:

```c
void sigma_xoshiro256ss_seed(sigma_xoshiro256ss_t *, uint64_t seed);
uint64_t sigma_rng_next_u64(sigma_state_t *);
double sigma_rng_next_normal(sigma_state_t *);
void sigma_rng_fill_normal(sigma_state_t *, double *out, size_t n);
```

The normal generator uses Marsaglia's polar transform. Independent ideal uniform
U,V on [0,1) give u=2U−1, v=2V−1, and a=u²+v². Reject until 0<a<1, then return

```math
(Z_1,Z_2)=(u,v)\sqrt{-2\log(a)/a}.
```

The second draw is cached. PRISM fills a buffer by repeated
`sigma_rng_next_normal` calls and saves/restores the cache. Thus changing buffer
sizes preserves that scalar sequence, including after odd-length calls. Finite
uniform grids and finite-state pseudorandom generation approximate ideal
Gaussian inputs; moments and CDF checks do not prove exact Gaussianity.

The pinned upstream buffered normal implementation does not first consume a
cached spare when it enters its two-at-a-time loop. Mixing scalar and buffered
calls can therefore produce a different sequence from repeated scalar calls.
A test reproduces this behavior independently of the adapter. That observation
establishes a sequence-contract limitation, rather than by itself demonstrating
an incorrect marginal distribution.

Upstream SIMD uniform fills interleave jumped streams and advance their base
state by jumps at each fill call. They do not have the scalar fill's sequence or
chunking policy. The audit checks this when the compiler selects a SIMD table
and records when no SIMD path was exercised. The local aarch64 audit exercised
NEON; AVX2 execution remains unvalidated locally. PRISM deliberately selects the
scalar table and scalar normal routine, irrespective of available uniform-fill
SIMD support. Throughput numbers for this adapter must not be labeled SIMD
normal performance.

## Tests, measurements, and limits

The provider tests check integer transitions, capacity and overlap handling,
repeatability, cache/chunk invariance, callback failure propagation, and
transactional covariance output. Normal checks use 202,000 draws and theoretical
six-standard-error smoke thresholds for raw moments of orders one through four
and nine fixed normal-CDF thresholds. A covariance integration smoke check and
logical-sample concurrency checks exercise the actual adapter. These are finite
internal checks rather than an exhaustive PRNG assessment.

The shared covariance implementation additionally checks lengths 257,511,513
at H=0.1,0.3,0.5,0.7,0.9 on the FFTW path. All innovation basis vectors recover
the numerical square root, whose row inner products are compared with dense
Toeplitz entries on both boundary rows and every diagonal. This deterministic
check has a separate absolute tolerance from the 4,096-draw lag-covariance
experiment. Raw lag estimates include approximate pointwise 95% intervals
using the known Gaussian entry variance; they are not simultaneous bands.
The dependency-free path additionally checks a dense n=513, H=0.7 factor column.
It does not run that entire large-length FFTW experiment.

`bench_innovations` measures checked normal-buffer fill and a full pipeline:
correlated-noise sampling, additive perturbation, and deterministic heat-bank
decomposition including diagnostics. Setup and three warmups are excluded from
steady-state timings. Normal fill requests n innovations; circulant sampling
requests 2(n−1), so they are distinct timing workloads. PCG and SIGMA use different
sequences and normal transforms. Timing alone does not establish statistical
quality or decomposition improvement. The archived observations are uncontrolled
smoke measurements, not a general relative-speed conclusion or a performance
study across independent sessions. The historical RNG, other distributions,
mergeable statistics, and SIMD normal implementations are not benchmarked here.

A separate optional CI configuration checks out the pinned SIGMA revision and
runs its integration tests. Local results do not establish that remote CI has
passed. Other optional backend adapters, covariance recalibration after nonlinear
marginal transforms, antithetic/common-random-number comparisons, and auxiliary
operators remain later work. Deterministic decomposition remains the default.
