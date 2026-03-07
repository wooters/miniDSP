# miniDSP

A small C library for audio DSP: signal measurement, FFT spectrum analysis, biquad filtering, GCC-PHAT delay estimation, and audio I/O.

## Core Engineering Principles

1. **Clarity over cleverness** — Write code that's maintainable, not impressive
2. **Explicit over implicit** — No magic. Make behavior obvious
3. **Composition over inheritance** — Small units that combine
4. **Fail fast, fail loud** — Surface errors at the source
5. **Delete code** — Less code = fewer bugs. Question every addition
6. **Verify, don't assume** — Run it. Test it. Prove it works

## Memory management

- **FFT plan caching** — `src/minidsp_spectrum.c` and `src/minidsp_gcc.c` keep static FFTW plans and buffers, allocated on first use and reallocated only when signal length changes. Call `MD_shutdown()` to free them.
- **Mel filterbank caching** — Filterbank matrices are cached by `(N, sample_rate, num_mels, min_freq_hz, max_freq_hz)` and released from `MD_shutdown()`.

## API design contract

- **Assertions vs sentinel returns** — Use `assert()` for structural misuse; sentinel returns (e.g. `0.0`) only for valid-but-unresolved runtime outcomes.
- **Module boundaries** — FFT-dependent APIs belong in `src/minidsp_spectrum.c`; time-domain/stateless analysis belongs in `src/minidsp_core.c`.
- **MFCC contract** — HTK mel mapping, one-sided PSD mel energies, natural-log floor at `1e-12`, DCT-II with `sqrt(1/M)` for `C0` and `sqrt(2/M)` for higher coefficients, `C0` returned at index 0.
- **Spectrum frequency bounds** — Assert structural misuse only; clamp runtime frequency ranges to `[0, Nyquist]`; return finite outputs for valid-but-empty clamped ranges.

## Documentation conventions

- **"Reading the formula in C" sections** — Use direct loops/arithmetic (not library wrappers) with explicit variable mapping comments (e.g. `i -> n`, `out[i] -> w[n]`).
- **Comparative visuals** — Keep plot settings fixed across variants (same tap length, FFT visualization length, dB floor/range) for meaningful comparison.
- **Guide assets** — When adding generated HTML plots/audio, update `Doxyfile` `HTML_EXTRA_FILES` and ensure every embedded iframe source is listed there.

## Build quirks

- **C17 standard** — The codebase targets `-std=c17 -Wall -Wextra -pedantic`. Use `NULL` (not `nullptr`), `#include <stdbool.h>` for `bool`, and `__attribute__((deprecated))` (not `[[deprecated]]`).
- **Legacy demo programs** in `tests/` (`gcc_phat_test`, `testliveio`, etc.) require gnuplot_i and/or PortAudio at link time — unlike the main test suite which only needs FFTW3, math, and libsndfile.
