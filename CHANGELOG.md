# Changelog

## v0.2.0-dev — unreleased

- Caller-owned PCG RNG with explicit seeds/streams, 53-bit uniforms, and cached
  Box–Muller normals; PCG attribution and Apache-2.0 license retained.
- Swappable caller-owned normal-source callbacks and a deterministic covariance
  transform API; the PCG convenience path preserves its previous draw sequence.
- Optional, header-hash-pinned SIGMA scalar normal adapter with explicit seed and
  logical-sample policy, checked chunk invariance, and no upstream types or flags
  exposed to consumers. Buffered/SIMD upstream sequence limitations are audited.
- Covariance square-root and lag-interval checks at lengths 257,511,513, plus
  serial/1/2/4/7-thread checks and normal-fill/end-to-end smoke benchmarks.
- Stable fGn lag covariance with explicit Hurst parameter and increment variance.
- Dense Cholesky, direct-DFT circulant reference, and optional FFTW circulant
  samplers with immutable plans, private workspaces, checked capacities, no
  PRISM sampling-time allocations, and mandatory factorization diagnostics.
- Explicit embedding rejection and opt-in, reported small-eigenvalue correction;
  no hidden regularization, automatic fallback, or Toeplitz/DCT substitution.
- Optional fixed-q ensemble experiment with 20 seeded replicates, zero/white/fGn
  controls, E=1,4,16,64,256, fixed amplitude, analytic moments, raw CSV, and a
  descriptive Monte Carlo uncertainty report.
- Covariance/oracle, allocation failure, and concurrent logical-stream tests;
  installed-consumer checks include the new headers and symbols.
- CMake and shared-library minor series moves to 0.2. Numeric package metadata is
  0.2.0; the public version string is 0.2.0-dev until a release is prepared.

Deterministic decomposition remains the default. This work checks numerical
covariance and fixed-linear identities; it does not demonstrate denoising,
source separation, inference calibration, or an ensemble benefit over the
noise-free deterministic operator. Other backend adapters, auxiliary operators,
and a production ensemble API remain outside this development scope.


## v0.1 — 0.1.0 (2026-10-21)

### Implemented

- C11 compiled library with orthonormal DCT analysis/synthesis contracts.
- Fixed nonoverlapping bands and telescoping heat partitions with an explicit
  residual, caller-specified parameters, and reconstruction diagnostics.
- Reference and optional FFTW backends, allocator callbacks, immutable plans,
  workspace-local execution state, and checked buffer capacities/overlap.
- Static/shared CMake packaging, installed-consumer checks, examples, a parameterized
  POSIX benchmark, and GCC/Clang CI configurations.
- Internal cosine-matrix oracle, deterministic property and operator tests,
  allocation/backend failure cases, and concurrent-workspace checks.

### Evidence and limits

The README specifies the mathematics, thesis, and null hypotheses. Numerical
results are local observations on the configurations in `docs/VALIDATION.md`.
The committed benchmark CSV files are preliminary smoke observations with unequal
sizes and uncontrolled hardware conditions. They establish no speed comparison.
Remote CI has not yet supplied release evidence.

This release provides deterministic global spectral filtering. It makes no
novelty, physical source recovery, application improvement, probabilistic
inference, or production-readiness claim. Adaptive/local methods, stochastic
ensembles, covariance sampling, auxiliary operators, and state-space models
remain future research and are not included in v0.1.
