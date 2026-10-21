# Validation record

## v0.2.0 development: innovation providers and SIGMA

The next development iteration added a caller-owned normal-source boundary,
a deterministic covariance-transform entry point, and the optional pinned
SIGMA scalar adapter. Local checks used the same Linux/aarch64 Jetson host,
Clang 22.1.3, GCC 11.4.0, FFTW 3.3.8, and SIGMA revision
`6e9a1ef526365474738d80d4eb71fd945a13abd6`.

| Configuration | Result |
| --- | --- |
| Clang Release static, FFTW and SIGMA enabled, experiments enabled | 19/19 checks passed |
| Clang Release shared, FFTW and SIGMA enabled, experiments enabled | 19/19 checks passed |
| GCC Release reference-only, SIGMA disabled, experiments enabled | 14/14 checks passed |
| GCC Debug, ASan/UBSan with leak detection, FFTW and SIGMA enabled | 19/19 checks passed; no sanitizer findings |
| GCC Release reference-only with SIGMA enabled | Provider, concurrency, SIGMA ensemble, and installed-consumer checks passed (4 selected checks) |
| Changed dependency header | Configuration rejected the altered pinned RNG header as expected |
| Installed C++11 consumer | New provider/SIGMA public headers compiled and linked with g++; no upstream public-header dependency |
| Clang static analysis | No diagnostics for innovation, SIGMA adapter, or modified covariance implementation |
| Markdown math | All 183 expressions passed source checks and MathJax 3.2.2 rendering |

A later provider test added finite-input overflow and verified that covariance
output stays unchanged. Its targeted static, shared, reference-only, and
sanitizer checks passed. The implementation itself was unchanged by that added
test. Sanitizer execution ran outside the sandbox for LeakSanitizer process
inspection. Full-suite counts above precede that targeted test extension.

The adapter consumes only the hash-checked upstream RNG/core headers and license;
its build imports no upstream fast-math or ISA options. The audit reproduced a
cached scalar/buffered normal sequence disagreement and NEON uniform-fill
sequence/chunking disagreement. The adapter deliberately uses scalar cached
normal calls; those observations do not by themselves establish a marginal
normal-distribution defect. AVX2 execution has not been validated locally.

Provider tests verify checked capacities and overlaps, callback error/nonfinite
propagation, unchanged covariance output on failure, explicit innovation counts,
and equivalence of the original PCG convenience path with the new source path.
SIGMA checks compare integer transitions against a separate oracle, preserve
cache behavior across odd chunks, and exercise 202,000 normals with moment and
fixed-grid CDF smoke checks. The maximum CDF standardized defect was 2.53024.
Serial and 1,2,4,7-thread sampling agree bitwise within the tested builds when
stream IDs are assigned to logical samples. These are finite internal checks,
not independence, cryptographic, or all-platform guarantees.

All innovation basis vectors recover each FFTW covariance square root at
n=257,511,513 and H=0.1,0.3,0.5,0.7,0.9. Checked boundary-row and diagonal inner
products differed from dense Toeplitz targets by at most 6.94e-16. The separate
4,096-draw lag experiment produced 6,405 estimates, with maximum standardized
defect 3.58148; 314 approximate pointwise 95% intervals excluded their targets.
Dependent lag estimates do not establish simultaneous coverage or interval
calibration. The large Monte Carlo draws use PCG, while SIGMA has its own
provider, integration, concurrency, and fixed-linear experiment checks.
A dependency-free dense n=513, H=0.7 factor-column check also passed.

[Raw provider/SIGMA results](../results/2026-10-21-providers-sigma/README.md)
include the upstream audit, all lag intervals, 400 SIGMA ensemble observations,
source hashes, and matched-size normal-fill/full-pipeline smoke timings.
The observed SIGMA ensemble mean/prediction ratios were 0.948–1.069 for nonzero
controls. No benefit over deterministic Tx is inferred. Timing is one
uncontrolled session, not a comparative performance conclusion. Other
distributions, statistics adapters, antithetic/common-random-number comparisons,
other backend adapters, and auxiliary operators remain unfinished. A pinned
SIGMA CI configuration is present; no remote CI success is claimed here.

## Earlier v0.2.0 development: covariance and ensemble

Local checks on 2026-10-21 used the Linux/aarch64 Jetson host described below,
Clang 22.1.3, GCC 11.4.0, and FFTW 3.3.8. This records a development snapshot,
not a published release or remote CI result.

| Configuration | Result |
| --- | --- |
| Clang Release, static, FFTW, experiments enabled | 14/14 CTest checks passed |
| Fresh Clang Release, shared, FFTW, experiments enabled | 14/14 CTest checks passed, including installed new-module consumer |
| GCC Release, reference-only, experiments enabled | 12/12 checks passed; added ownership/failure cases also passed after the full run |
| GCC Debug, ASan/UBSan and leak detection, FFTW, experiments enabled | 14/14 checks passed across the suite and the final clean-source rerun; no sanitizer findings |
| C++11 stochastic-header/link smoke check | Complete noise example compiled and ran with g++ |
| New implementation static analysis | Clang analyzer reported no diagnostics for RNG and noise units |
| Mathematics source/render checks | 172 expressions across published Markdown passed the source checker and MathJax 3.2.2 rendering |

The sanitizer suite ran outside the sandbox because LeakSanitizer requires
process inspection. Its clean-source check initially copied the developing
results README before the referenced source-hash artifact existed; that
packaging/documentation check passed after the artifact was completed. The
instrumented execution checks passed without sanitizer findings. A reused shared
build cache after a compiler change requested a missing CMake relink artifact;
the fresh shared build and installation passed. Neither transient failure is
excluded from this record or treated as a scientific result.

New tests check the published PCG integer vector, cached/chunked normal draws,
fGn covariance against direct power differences at small lags, sign preservation
near H=1/2, and exact white/zero covariance controls. An independent full
cosine-matrix circulant square root checks normalization and outputs at lengths
1,2,3,8,17. A positive definite Toeplitz target with an indefinite mirror
embedding is accepted by the dense method and rejected by both circulant
methods. Explicit small-negative rounding is separately checked and reported.

Known-zero-mean empirical second moments use 8,192 draws at each of
H=0.1,0.3,0.5,0.7,0.9, with dense n=32 and FFTW circulant n=128. Six selected
entries are compared with Gaussian theoretical standard errors, and all entries
contribute to the Frobenius defect. Maximum selected-entry standardized defect
was 2.0458; full-matrix relative defects ranged from 0.0243 to 0.1247.
This is finite-sample internal evidence, not an independence proof or a uniform
sampler error bound. Raw diagnostics, protocol, and limitations are in the
[covariance/ensemble results](../results/2026-10-21-covariance-ensemble/README.md).

Allocation-failure injection checks cleanup for each method; a locked allocator
rejects explicit execution-time allocations. Tests also check capacities,
sentinels, copied configuration, invalid RNGs, output/state overlap, workspace
ownership, and serial versus four-thread logical-stream scheduling. FFTW
planner creation/destruction is serialized in those tests. This is not a formal
race analysis or allocator interposition inside external libraries.

Both ensemble drivers ran the fixed protocol at 20 seeded replicates per control
with E=1,4,16,64,256. The maximum paired linear-identity defect was 1.67e-15.
Squared-error means passed the six-theoretical-standard-error smoke threshold.
The white E=1 approximate 95% t interval excluded the predicted mean in both
runs; intervals are descriptive and small-sample coverage is not guaranteed.
Deterministic Tx already has zero perturbation error. These checks establish no
held-out application result or ensemble benefit over that baseline.

Reproduce the development checks by adding `-DPRISM_BUILD_EXPERIMENTS=ON` to the
build commands below. The experiment commands, numeric protocol, source hashes,
and raw observations are linked above. Normal samples are not promised bitwise
identical across platforms. Remote CI, larger sampler sizes, and other platforms
remain unestablished by this local record.

## v0.1 validation archive

The initial v0.1 checks on 2026-10-21 used an NVIDIA Jetson Orin Nano Engineering Reference
Developer Kit Super, Linux 5.15.148-tegra aarch64, FFTW 3.3.8, Clang 22.1.3, and
GCC 11.4.0. All project targets use C11 and strict warnings, including conversion
warnings, treated as errors. Release tests retain always-active CHECK macros.

| Configuration | Result |
| --- | --- |
| Clean Clang Release, static library, FFTW enabled | 10/10 CTest checks passed, including the documentation gate |
| Clean Clang Release, shared library, FFTW enabled | 9/9 CTest checks passed |
| GCC Release, reference-only | 8/8 CTest checks passed |
| GCC Debug, ASan/UBSan and leak detection, FFTW enabled | 9/9 CTest checks passed; no sanitizer findings |
| Installed CMake consumer | Static/shared and reference-only packages passed, including the README C example and incompatible-version rejection |
| Subproject consumer | Parent C standard preserved; optional targets defaulted off |
| Clean-source snapshot | Builds and selected checks passed using only Git-includable files |
| C++11 public-header/link smoke check | Passed with g++ |
| Clang static analyzer, all three implementation units | No diagnostics |
| Retired-name content/filename audit | Zero matches in project files, local configuration, and fresh build artifacts |

The documentation gate was added after the shared/reference/sanitizer runs
recorded above. With Python 3 available, it also runs inside the source snapshot,
where local planning files are absent. The reported CI failure was reproduced
from the committed tree and the corrected checker passed that same file set.

The local sandbox prevents LeakSanitizer process inspection. Sanitizer tests
were rerun outside it with `ASAN_OPTIONS=detect_leaks=1:halt_on_error=1` and
`UBSAN_OPTIONS=halt_on_error=1`. Sandbox inspection failures are not numerical or
memory-correctness findings. CI runs the GCC/Clang, FFTW/reference-only matrix
and an additional sanitizer and shared-library configurations. Remote CI has not yet
provided release evidence.

### Numerical coverage

Both backends and both partitions passed reconstruction checks at lengths
128,257,511,512,513,1024,4096,16384. The largest observed relative reconstruction
error was approximately 2.00e-13 for the reference backend and 5.10e-16 for FFTW,
below the 1e-10 test target. Partition defect was below 1e-14.

Transform tests compare DC, impulses, sinusoids, and every basis vector at a
small length against an independently written direct cosine-matrix oracle.
Component outputs, rather than just their sum, are compared with that oracle.
Zero, constant, boundary impulse, highest-frequency basis, and scales 1e-150 and
1e150 are exercised at lengths through 1024; large lengths use a mixed trend and
oscillation. The reference backend uses periodically re-anchored trigonometric
rotation and compensated sums; the oracle uses direct cosine evaluations.

An additional 100 seeded cases per backend cover varied lengths, spacing,
linearity, DC shifts, reconstruction, and energy including cross terms. Other
tests cover checked capacities and sentinel preservation, partial overlaps,
configuration ownership, size overflow, nonfinite input, allocation failures,
backend errors, numerical overflow, repeated execution, and simultaneous use of
one plan with distinct workspaces. A tracked allocator rejects execution-time
allocation and checks that all caller-allocator allocations are released.

The mathematical-property checks additionally verify heat semigroup products,
physical-unit rescaling, positive close-time bands, fixed-band idempotence, and
per-component/stacked L2 bounds against the README definitions. The installed
consumer compiles and runs the complete README C block, accepts the 0.1 package,
and rejects a 0.2 request. Subproject and clean-source checks exercise packaging
without relying on ignored local files. These are internal checks, not external
peer review or formal verification.

The local-link checker and manual CLI smoke checks cover malformed backend names,
negative/overflowing dimensions, zero repetitions, and matched-size successful
benchmark calls. The archived CSV files remain unchanged; their unconditional
preliminary planner label is explained in the results note.

### Markdown rendering check

After correcting the initial GitHub formatting, all 113 README equations were
submitted to GitHub's Markdown preview endpoint. The TeX inside its math-renderer
elements matched the protected source exactly, and those expressions rendered
without errors using MathJax 3.2.2 in a temporary local validation environment.
That renderer is not a PRISM build or CI dependency. This was a preview/source
and renderer check, rather than an end-to-end browser test.

The CI documentation checker now flags the known unsupported operator macro,
unprotected math delimiters, and unbalanced TeX braces. It rejects the original
README that exposed the problem. Other published Markdown files contained no
math expressions. The source checker is not a complete TeX parser or proof of
compatibility with every future GitHub renderer.

### Reproduction

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel 3
ctest --test-dir build -L prism --output-on-failure --parallel 3
cmake -S . -B build-reference -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_BUILD_TYPE=Release -DPRISM_ENABLE_FFTW=OFF -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build-reference --parallel 3
ctest --test-dir build-reference -L prism --output-on-failure --parallel 3
cmake -S . -B build-sanitize -DCMAKE_C_COMPILER=gcc \
    -DCMAKE_BUILD_TYPE=Debug -DPRISM_SANITIZE=ON -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build-sanitize --parallel 3
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    ctest --test-dir build-sanitize -L prism --output-on-failure --parallel 3
```

[Local benchmark outputs](../results/2026-10-21-jetson-orin-nano/README.md) are one
uncontrolled smoke session, not a comparative performance result. Cross-platform
and application-level evidence remains future work. Historical failure
reproductions are described in [the external baseline audit](BASELINE-AUDIT.md).
