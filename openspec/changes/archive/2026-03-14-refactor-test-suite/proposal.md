## Why

`tests/test_minidsp.c` has grown to 5,757 lines and 277 tests in a single file. It's hard to navigate, slow to recompile when touching any test, and will only get worse as more API functions are added. Splitting it into per-module test files — mirroring the existing `src/minidsp_*.c` layout — will make tests easier to find, faster to compile incrementally, and simpler to extend.

## What Changes

- Split `test_minidsp.c` into ~12 focused test files, one per logical module (core, generators, spectrum, gcc, biquad, fir, effects, steg, dtmf, spectext, fileio, resample).
- Extract the shared test infrastructure (counters, `approx_equal`, `RUN_TEST` macro, `delay_signal` helper) into a shared header (`test_helpers.h`).
- Each test file exposes a single `void run_xxx_tests(void)` function that prints its section headers and calls `RUN_TEST` for its group.
- A slim `test_minidsp.c` remains as the driver: includes all test files, calls each `run_xxx_tests()`, prints the summary, calls `MD_shutdown()`.
- Update `tests/Makefile` to compile all test `.c` files and link them together.
- No test logic changes — all 277 tests remain identical in behavior.

## Capabilities

### New Capabilities
- `split-test-files`: Splitting the monolithic test file into per-module files with a shared harness and driver.

### Modified Capabilities
_(none — no spec-level behavior changes)_

## Impact

- **Files changed**: `tests/test_minidsp.c` (gutted to driver), `tests/Makefile` (multi-object build).
- **Files added**: `tests/test_helpers.h`, `tests/test_core.c`, `tests/test_generators.c`, `tests/test_spectrum.c`, `tests/test_gcc.c`, `tests/test_biquad.c`, `tests/test_fir.c`, `tests/test_effects.c`, `tests/test_steg.c`, `tests/test_dtmf.c`, `tests/test_spectext.c`, `tests/test_fileio.c`, `tests/test_resample.c`.
- **CI/containers**: `make test` still works identically — just compiles more objects.
- **No API or library changes.**
