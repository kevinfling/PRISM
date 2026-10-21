# Contributing to PRISM

PRISM v0.1 implements the deterministic operators specified in the README.
Keep changes small enough to review against that specification. The public
API is `include/prism.h`; implementation is in `src/`. Use the repository's
`.clang-format` style and preserve copyright and license notices.

## Checks

Build with strict warnings, run the complete suite, then use a separate build
for reference-only and sanitizer checks:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPRISM_WARNINGS_AS_ERRORS=ON
cmake --build build --parallel
ctest --test-dir build -L prism --output-on-failure
cmake -S . -B build-reference -DPRISM_ENABLE_FFTW=OFF -DCMAKE_BUILD_TYPE=Release
cmake --build build-reference --parallel
ctest --test-dir build-reference -L prism --output-on-failure
cmake -S . -B build-sanitize -DCMAKE_BUILD_TYPE=Debug -DPRISM_SANITIZE=ON
cmake --build build-sanitize --parallel
ASAN_OPTIONS=detect_leaks=1:halt_on_error=1 UBSAN_OPTIONS=halt_on_error=1 \
    ctest --test-dir build-sanitize -L prism --output-on-failure
```

The installed-consumer test checks the exported CMake package. The clean-source
check exercises only files Git would include in the initial commit, so ignored
local files cannot mask a missing dependency. When Python 3 is available, CTest
also runs the documentation checker in both the working tree and source snapshot.
The CI configuration also includes shared-library and
subproject integration checks. Run `python3 tools/check_docs.py` to check local
Markdown links and anchors (Python 3 standard library only). Test assertions must remain active under NDEBUG.

Add tests when behavior, numerical definitions, or ownership contracts change.
Use analytic cases or a separately implemented oracle rather than duplicating
the implementation. A reconstruction-only check can miss incorrect individual
bands. Failure-path tests should check capacities, partial overlap, cleanup,
and whether result metadata is usable after an error.

## Markdown mathematics

Use fenced `math` blocks for displayed equations and dollar/backtick delimiters
for inline equations, as described in [GitHub's math documentation](https://docs.github.com/en/get-started/writing-on-github/working-with-advanced-formatting/writing-mathematical-expressions).
These protect TeX norm bars, braces, and spacing commands from Markdown escaping.
Use `\mathop{\mathrm{diag}}`-style operators instead of `\operatorname`, which
GitHub rejects. Prefer `\lt` and `\gt` for comparisons, and keep a TeX row break
on the same source line as the following row to avoid an added Markdown line break.

`python3 tools/check_docs.py` catches the known macro and delimiter pitfalls as
well as local link errors. It is not a complete TeX parser; review the rendered
GitHub preview when equations change. Source math can be valid TeX yet be damaged
by the Markdown stage before it reaches the renderer.

## Mathematical and empirical evidence

Document sample positions, transform normalization, units, boundaries, and the
residual before introducing a new operator. Separate exact-arithmetic identities
from numerical measurements. A separately written test oracle is an internal
cross-check, not independent external validation. Future stochastic code needs
explicit seeds, covariance definitions, and controls; do not add hidden RNG state.

For a research claim, define the task, null hypothesis, primary metric, minimum
useful effect, matched resource budget, held-out split, and analysis before final
evaluation. Archive raw observations, commands, software versions, and machine
conditions. Report uncertainty and negative results. Do not infer separation,
independence, causality, uncertainty calibration, or scientific novelty from
reconstruction accuracy alone.

Benchmark both backends at the same length and repetition count. Distinguish
setup, first execution, warm execution, and complete application cost. One
uncontrolled smoke session is not a comparative performance result.

## Compatibility and review

Public changes during the 0.1 series should preserve the documented API where
practical. Later pre-1.0 minor releases may break compatibility; record those
changes in the changelog. Keep optional backend details out of public headers.
Plans must outlive workspaces; each concurrent execution uses a distinct
workspace. Borrowed callback descriptors and allocator state must obey the
lifetimes in `docs/API.md`.

For a proposed change, explain the user-visible behavior, the mathematical or
engineering reason, relevant test evidence, and any limits. Update README/API
contracts together with code. Proposed research does not constitute a list of
shipped capabilities.
