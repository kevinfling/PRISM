# Preliminary local benchmark smoke results

These files were recorded with the preliminary driver before release polish;
they are retained unchanged. Commands from the project root, after a Clang 22.1.3
Release build:

```sh
./build/benchmarks/bench_prism > results/2026-10-21-jetson-orin-nano/fftw.csv
./build/benchmarks/bench_prism reference > results/2026-10-21-jetson-orin-nano/reference.csv
```

Machine: NVIDIA Jetson Orin Nano Engineering Reference Developer Kit Super,
Linux 5.15.148-tegra aarch64. Compiler options: C11, Release -O3 -DNDEBUG,
-Wall -Wextra -Wpedantic -Wconversion -Werror; no forced native CPU/fast-math.
The FFTW row uses FFTW 3.3.8, ESTIMATE|UNALIGNED, single-threaded.
The reference row uses the quadratic cosine implementation; its printed planner
metadata was unconditional in that preliminary driver and does not describe an
FFTW execution. The current driver reports backend-specific metadata. PRISM buffers use malloc/free;
FFTW planning accepts arbitrary ordinary-double alignment. No RNG is used.

Each file records plan/workspace setup, first execution, and the mean of 20
executions after three warm-up calls. Timings include decomposition validation,
transforms, aggregation, and reconstruction diagnostics. They exclude signal
creation, caller output allocation, plan destruction, and console output.

This is one smoke session with uncontrolled CPU frequency, background load,
and power settings. Backend rows use different input lengths; do not interpret
them as a speed comparison or a release performance claim. Three independent
sessions with matched sizes, planner effort, and controlled hardware conditions
are required for a release comparison. The reference implementation is quadratic.

For new matched-size observations, use explicit arguments, for example
`bench_prism fftw 512 20` and `bench_prism reference 512 20`, and save the new
outputs separately. The current driver rejects an unavailable requested backend
and malformed dimensions. A mean warm time does not describe a latency tail.
