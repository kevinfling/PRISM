# v0.1 validation record

Local checks on 2026-10-21 used an NVIDIA Jetson Orin Nano Engineering Reference
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

## Numerical coverage

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

## Markdown rendering check

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

## Reproduction

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
