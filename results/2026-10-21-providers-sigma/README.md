# v0.2.0 development: innovation providers and SIGMA

Local observations from 2026-10-21 on the Linux/aarch64 NVIDIA Jetson Orin Nano
Engineering Reference Developer Kit Super, Clang 22.1.3 Release, FFTW 3.3.8,
and SIGMA revision `6e9a1ef526365474738d80d4eb71fd945a13abd6`.
This uncommitted development snapshot is based on v0.1 commit
`9fa3bb72d7555baa85491d55c84a4275d7e09d1b`. The [source hashes](source-sha256.txt)
identify this iteration; the [earlier covariance/ensemble session](../2026-10-21-covariance-ensemble/README.md)
retains its original observations and hashes. No release, external validation,
application benefit, or general PRNG-quality claim follows from this session.

See the [adapter audit and protocol](../../docs/SIGMA.md) and
[provider/covariance contracts](../../docs/STOCHASTIC.md#innovation-providers-and-deterministic-covariance-transforms).

| Artifact | Observation |
| --- | --- |
| [normal-checks.txt](normal-checks.txt) | 202,000 SIGMA scalar cached normals; raw moments and nine normal-CDF threshold checks |
| [upstream-audit.txt](upstream-audit.txt) | Cached scalar/buffered normal sequence mismatch; NEON-selected uniform sequence and chunking disagreement |
| [covariance-checks.txt](covariance-checks.txt) | Deterministic square-root checks and lag-covariance summaries for H=0.1,0.3,0.5,0.7,0.9 at lengths 257,511,513 |
| [lag-covariance.csv](lag-covariance.csv) | All 6,405 first-row lag estimates, targets, theoretical standard errors, and approximate pointwise 95% intervals; 4,096 draws per case |
| [ensemble-sigma.csv](ensemble-sigma.csv), [ensemble-sigma.txt](ensemble-sigma.txt) | All 400 replicate/count observations and fixed-linear protocol metadata using the SIGMA source |
| [summary-sigma.csv](summary-sigma.csv) | Descriptive Monte Carlo means, empirical standard errors, approximate t intervals, and theoretical predictions |

The normal sample mean was −0.0028823; raw second and fourth moments were 1.00464
and 3.04302. The maximum fixed-grid CDF defect was 2.53024 theoretical standard
errors. Moment and CDF checks passed the six-standard-error smoke thresholds.
This is neither a comprehensive distribution test nor an independence proof.

For FFTW covariance plans, all innovation basis vectors recover a numerical
square root. Checked boundary-row and diagonal covariances differed from their
Toeplitz targets by at most 6.94e-16. The largest Monte Carlo lag defect was
3.58148 theoretical standard errors. Of 6,405 approximate pointwise intervals,
314 excluded their targets. The lag estimates are dependent within each case;
this count does not establish interval calibration or simultaneous coverage.
No eigenvalue rounding was used. A dense n=513, H=0.7 first-factor-column check
also passed. Large-size Monte Carlo draws here use the PCG source, while SIGMA
is exercised by the separate distribution, concurrency, integration, and
ensemble checks; do not attribute every result to SIGMA.

The SIGMA ensemble mean/prediction ratios ranged from 0.948 to 1.069 over
nonzero controls/counts. No theoretical mean lay outside its descriptive
approximate t interval in this particular run. Those finite observations do
not prove coverage or alter the earlier PCG interval exception. Deterministic
Tx already has zero perturbation error; these runs establish no ensemble or
denoising advantage over that baseline.

Normal-fill and full-pipeline smoke timings use the same sizes and 100 calls
per provider, with setup and three warmups recorded separately:

| n | PCG | SIGMA |
| --- | --- | --- |
| 2 | [raw](benchmark-pcg-n2.csv) | [raw](benchmark-sigma-n2.csv) |
| 3 | [raw](benchmark-pcg-n3.csv) | [raw](benchmark-sigma-n3.csv) |
| 257 | [raw](benchmark-pcg-n257.csv) | [raw](benchmark-sigma-n257.csv) |
| 1024 | [raw](benchmark-pcg-n1024.csv) | [raw](benchmark-sigma-n1024.csv) |

These are consecutive calls in one uncontrolled session, not independent
performance sessions or a general speed comparison. Clock/power state and host
load were not controlled. Normal fills include checked public-API overhead;
the full pipeline includes noise sampling, perturbation, decomposition, and
reconstruction diagnostics. The adapter uses scalar cached normals, not SIMD
normal fills. Different providers use different sequences and normal transforms.
No historical RNG, probability distribution suite, or mergeable online-statistics
adapter was measured.

Reproduce with the [pinned build commands](../../docs/SIGMA.md#dependency-boundary-and-pin), then:

```sh
./build-sigma/tests/test_innovation > normal-checks.txt
./build-sigma/tests/test_sigma_upstream > upstream-audit.txt
./build-sigma/tests/test_noise_extended lag-covariance.csv > covariance-checks.txt
./build-sigma/experiments/prism_ensemble_experiment fftw sigma > ensemble-sigma.csv 2> ensemble-sigma.txt
python3 tools/summarize_ensemble.py ensemble-sigma.csv > summary-sigma.csv
./build-sigma/benchmarks/bench_innovations pcg 257 100 > benchmark-pcg-n257.csv
./build-sigma/benchmarks/bench_innovations sigma 257 100 > benchmark-sigma-n257.csv
```

Repeat benchmark calls at n=2,3,1024 for the other rows. Normal values and timings
are not promised bitwise identical across builds or platforms. The local audit
exercised NEON uniform fills; AVX2 execution and remote CI are not established
by this record.
