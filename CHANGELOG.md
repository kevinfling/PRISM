# Changelog

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
