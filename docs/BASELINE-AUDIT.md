# External baseline audit

The unmodified project snapshot and SHA-256 manifest are stored outside this
repository at `/tmp/PRISM-baseline-2026-10-21/`. That temporary archive must be
copied to durable storage. It contains original license/attribution, source,
examples, tests, documentation, and local configuration.

This is a local historical audit. The archive is not distributed with PRISM,
and a fresh clone cannot reproduce these observations without obtaining it
separately. The PRISM correctness tests do not depend on the archive. The bundled
snapshot is `PRISM-baseline-2026-10-21.tar.gz`; its SHA-256 is recorded below so
that a separately retained copy can be identified.

```text
250b2e8c6cd1f265ac3499d825cd78cfb9d181777c0876eefc535bd64388f0ea
```

The archive's `diagnostics/` directory contains a standalone `legacy_behavior.c`,
pinned CHEAP header and license, commands, numerical vectors, and sanitizer logs.
The dependency is from commit `90e343a841f79893d2d7b7fb252e57110b2a0562`.
The harness includes the dependency first, then overrides the timestamp call with
seed 1234567 before including the original consumer header. None of the six
original headers was edited. Full-output repetitions matched on this machine.

| Finding | Evidence class | Result |
| --- | --- | --- |
| Missing output capacity | Reproduced with ASan | One-row output allocation overflowed during extraction. |
| Missing count-pointer validation | Reproduced with UBSan | NULL count pointer caused a null store. |
| NULL output shares a trial buffer | Source-observed and reproduced | Under identical seed/input, first residual sample was 2.139655883969148 instead of 2.1791139447115722 with full output. This establishes behavior divergence; neither vector defines a new correctness contract. |
| Noise scale depends on ensemble count | Source-observed | Noise is scaled by the inverse square root of count before averaging; statistical evaluation is deferred. |
| Exact target covariance | Unverified | A DCT-diagonal Gaussian field must not be assumed to have the desired Toeplitz fGn covariance. |
| Scalar-state model and oscillator placeholder | Source-observed | Neither module is exposed in v0.1; independent state-space reference tests are future work. |
| Missing log residual and fractional-definition ambiguity | Source-observed | These modules are archived and deferred pending explicit reconstruction/operator contracts. |

Reproduction from the external archive requires GCC, FFTW3, and libm:

```sh
cd /tmp/PRISM-baseline-2026-10-21/diagnostics
gcc -std=c11 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
    -I . -I ../include legacy_behavior.c -lfftw3 -lm -o legacy_behavior
python3 run.py
```

The runner saves successful numerical outputs separately from expected failures.
ASan/UBSan failures in these historical diagnostics are evidence, not PRISM test
expectations. Sanitizer runs need an environment that permits LeakSanitizer's
process inspection. Original SHA-256 hashes were verified after diagnostics.
