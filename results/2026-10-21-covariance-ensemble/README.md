# v0.2.0 development: covariance and fixed-linear ensemble checks

Local observations from 2026-10-21 on an NVIDIA Jetson Orin Nano Engineering
Reference Developer Kit Super, Linux/aarch64, Clang 22.1.3 Release, FFTW 3.3.8.
This is an uncommitted development snapshot based on v0.1 commit
`9fa3bb72d7555baa85491d55c84a4275d7e09d1b` when generated. The
[source hashes](source-sha256.txt) identify the implementation and protocol used.
These observations are internal numerical checks, not independent external
validation, held-out application results, or a v0.2.0 release claim.

The [fixed protocol and mathematical predictions](../../docs/STOCHASTIC.md)
were specified in source before these runs. Each backend run contains 20 seeded
replicates for each of four controls and five ensemble counts: 400 CSV rows.
Counts within a replicate use nested prefixes and are dependent.

| Artifact | Contents |
| --- | --- |
| [covariance.txt](covariance.txt) | 8,192 draws per H at H=0.1,0.3,0.5,0.7,0.9; dense n=32 and FFTW circulant n=128; selected-entry standardized defects and full-matrix Frobenius defects |
| [ensemble-fftw.csv](ensemble-fftw.csv) | All fixed-linear replicate/count observations using FFTW DCT and circulant noise |
| [ensemble-reference.csv](ensemble-reference.csv) | All observations using reference DCT and dense Cholesky noise |
| [summary-fftw.csv](summary-fftw.csv), [summary-reference.csv](summary-reference.csv) | Descriptive means, empirical Monte Carlo standard errors, approximate t intervals, and theoretical moments |
| [ensemble-fftw.txt](ensemble-fftw.txt), [ensemble-reference.txt](ensemble-reference.txt) | Version, fixed parameters, seeds, backend/planner policy, and factorization diagnostics |

The largest selected covariance entry defect was 2.0458 theoretical standard
errors; the smoke threshold is six. Full-matrix relative Frobenius defects were
0.0243–0.1247, with predicted RMS values 0.0241–0.1255 depending on H and size.
These quantify finite-sample covariance-estimation variability; they are not
deterministic factorization-error bounds. No eigenvalue rounding was used in these runs.

For the FFTW ensemble run, mean squared perturbation error divided by its
prediction ranged from 0.876 to 1.057 across nonzero controls and counts.
The white E=1 approximate 95% t interval did not contain the predicted mean,
although it passed the prespecified six-theoretical-standard-error
smoke check. The intervals are descriptive, unadjusted, and not guaranteed to
have nominal small-sample coverage. They must not be presented as uniformly
confirming the theory. Zero-control errors at large counts are roundoff from
repeated summation.

The expected error decreases as 1/E relative to the perturbed E=1 procedure.
The deterministic baseline already has zero perturbation error, so these runs
show no ensemble benefit over that baseline. They evaluate no denoising or
physical-source recovery task. Distinct PCG streams enable repeatable indexing;
these observations do not prove probabilistic independence or exact Gaussianity.

Component timings are uncontrolled smoke observations. Background validation
jobs were active during this session. Setup, one deterministic execution,
cumulative noise sampling, transforms, and aggregation are recorded separately;
oracle checks, diagnostics, and printing are excluded from those component
measurements. The noise algorithms differ between runs, so their times cannot
isolate DCT-backend speed. No comparative performance claim is supported.

Reproduce the protocol and reports:

```sh
cmake -S . -B build-polish -DCMAKE_C_COMPILER=/usr/bin/clang-22 \
    -DCMAKE_BUILD_TYPE=Release -DPRISM_WARNINGS_AS_ERRORS=ON \
    -DPRISM_BUILD_EXPERIMENTS=ON
cmake --build build-polish --parallel 3
./build-polish/tests/test_noise > covariance.txt
./build-polish/experiments/prism_ensemble_experiment fftw > ensemble-fftw.csv 2> ensemble-fftw.txt
./build-polish/experiments/prism_ensemble_experiment reference > ensemble-reference.csv 2> ensemble-reference.txt
python3 tools/summarize_ensemble.py ensemble-fftw.csv > summary-fftw.csv
python3 tools/summarize_ensemble.py ensemble-reference.csv > summary-reference.csv
```

Adjust the compiler path for other hosts. Cross-platform normal-sample bitwise
identity is not promised. The report generator uses only the Python standard
library; timing reruns will differ even on this host.
