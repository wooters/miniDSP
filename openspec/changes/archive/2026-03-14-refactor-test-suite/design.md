## Context

`tests/test_minidsp.c` is a single 5,757-line file containing all 277 tests for every module in the library. The library source is already cleanly split across 9 `src/minidsp_*.c` files plus `biquad.c` and `fileio.c`. The test file has no such structure — it's a flat sequence of static functions with a monolithic `main()`.

The test framework is home-grown: global counters, an `approx_equal` helper, and a `RUN_TEST` macro. There are no external test framework dependencies, and we want to keep it that way.

## Goals / Non-Goals

**Goals:**
- Split tests into per-module files that mirror `src/` layout for navigability.
- Enable incremental recompilation — changing one test file only recompiles that `.o`.
- Make it obvious where to add new tests for a given API function.
- Preserve the exact same 277 tests, output format, and exit behavior.

**Non-Goals:**
- Adopting an external test framework (CUnit, Check, Unity, etc.).
- Adding new tests or changing existing test logic.
- Changing the test output format or return codes.
- Parallel test execution or test isolation.

## Decisions

### 1. File-per-module split (not file-per-function or file-per-section)

**Decision**: One test file per source module, matching the `src/minidsp_*.c` naming.

**Rationale**: The source already has a clear module boundary (core, generators, spectrum, gcc, fir, resample, dtmf, spectext, steg). Mirroring this in tests creates a predictable mapping: "where do I test `MD_foo`?" → look at the same module as `src/`. A finer split (per-function) would create too many tiny files; a coarser split would still leave large files.

**Alternatives considered**:
- *Single file with `#include` splitting*: Fragile, misleading (looks like one TU but isn't), confusing for editors.
- *One file per test section header*: ~40 files, too granular for the project size.

### 2. Shared header for test infrastructure

**Decision**: `tests/test_helpers.h` contains the global counters, `approx_equal`, `RUN_TEST` macro, and the `delay_signal` helper. Each test file `#include`s it.

**Rationale**: The counters must be shared across all files for the summary to work. Making them `extern` in a header with a single definition in the driver keeps the framework minimal.

**Detail**: `test_helpers.h` declares the counters as `extern`. The driver (`test_minidsp.c`) defines them. The `RUN_TEST` macro and `approx_equal` stay as a macro and `static inline` respectively (no linkage issues).

### 3. Each test file exposes one registration function

**Decision**: Each `test_xxx.c` file contains its `static` test functions plus one non-static `void run_xxx_tests(void)` function that prints section headers and calls `RUN_TEST`.

**Rationale**: Keeps each file self-contained. The driver just calls `run_xxx_tests()` for each module. Test functions remain `static` (internal linkage), avoiding symbol collisions if two modules use similar helper names.

### 4. Driver file stays named `test_minidsp.c`

**Decision**: The existing `test_minidsp.c` becomes a slim driver that defines the global counters, includes headers for each test module's `run_*_tests()`, calls them, prints the summary, and calls `MD_shutdown()`.

**Rationale**: Preserves the existing `make test_minidsp` target and CI integration. No Makefile target renames needed.

### 5. Makefile uses multi-object link

**Decision**: Each `test_*.c` compiles to `test_*.o`. The `test_minidsp` binary links all `.o` files together.

**Rationale**: Standard C multi-file compilation. Incremental rebuilds only recompile changed `.o` files.

### File mapping

| Test file | Source module | Sections covered |
|---|---|---|
| `test_core.c` | `minidsp_core.c` | dot, energy, power, power_db, scale, fit_within_range, adjust_dblevel, entropy, rms, zcr, autocorrelation, peak_detect, pitch, mix |
| `test_effects.c` | `minidsp_core.c` (effects) | delay_echo, tremolo, comb_reverb |
| `test_generators.c` | `minidsp_generators.c` | sine, white_noise, impulse, chirp_linear, chirp_log, square, sawtooth, shepard_tone |
| `test_spectrum.c` | `minidsp_spectrum.c` | magnitude_spectrum, PSD, phase_spectrum, STFT, mel/MFCC |
| `test_fir.c` | `minidsp_fir.c` | convolution, FIR filter, moving average, FFT OLA, bessel_i0, sinc, Kaiser window, lowpass FIR design |
| `test_gcc.c` | `minidsp_gcc.c` | GCC-PHAT delay estimation |
| `test_biquad.c` | `biquad.c` | Biquad filter types |
| `test_windows.c` | `minidsp_core.c` (windows) | Hann, Hamming, Blackman, Rect windows |
| `test_dtmf.c` | `minidsp_dtmf.c` | DTMF generation and detection |
| `test_spectext.c` | `minidsp_spectext.c` | Spectrogram text |
| `test_steg.c` | `minidsp_steg.c` | Audio steganography |
| `test_fileio.c` | `fileio.c` | NPY, safetensors, WAV writers |
| `test_resample.c` | `minidsp_resample.c` | Resampling output length, resampling |

## Risks / Trade-offs

- **[Risk] Symbol collisions** → All test functions remain `static` within their file. Only the `run_*_tests()` functions have external linkage. Collision risk is minimal.
- **[Risk] Shared state via global counters** → The counters (`tests_run`, `tests_passed`, `tests_failed`) are `extern` globals defined in the driver. This is the simplest approach and matches the current design. Not a concern for a single-threaded test suite.
- **[Risk] Large diff** → This is a pure restructuring with no logic changes. Reviewing can be done by confirming: (a) every `RUN_TEST` from the original appears exactly once in the new files, (b) `make test` still passes with 277/277.
- **[Tradeoff] `test_core.c` is still large** → It covers ~14 sections because `minidsp_core.c` has many functions. This is acceptable since it mirrors the module boundary. If `minidsp_core.c` is ever split, the tests can follow.
