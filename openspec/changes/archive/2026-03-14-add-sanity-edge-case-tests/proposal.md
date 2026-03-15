## Why

The test suite has strong "functional correctness" tests — verifying results against known mathematical properties — but inconsistently covers trivial/degenerate inputs. Some functions (e.g., `MD_dot`, `MD_energy`, `MD_impulse`) have thorough sanity checks (zeros, single element, negative amplitude), while others lack obvious cases like zero-length-ish inputs, identity operations (e.g., delay=0, depth=0, width=1), or constant signals. Adding systematic sanity and edge-case tests catches regressions where a function might return NaN, inf, or garbage on trivial inputs that "should obviously work."

## What Changes

- Add **sanity tests** (trivial/obvious inputs with known outputs) to functions that currently lack them: zero signal, constant signal, single element, identity parameters.
- Add **boundary/edge-case tests** for degenerate but valid parameter values: zero amplitude, zero delay, width=1, feedback=0, same start/end frequency, single-frame STFT, length-1 kernels, same-rate resampling.
- All new tests follow the existing convention: `static int test_xxx(void)` returning 1=pass/0=fail, registered with `RUN_TEST()` in the appropriate section of `main()`.
- No changes to library code — this is purely additive test coverage.

## Capabilities

### New Capabilities
- `sanity-edge-case-tests`: Systematic sanity and edge-case test coverage for all public API functions in `tests/test_minidsp.c`

### Modified Capabilities

## Impact

- **`tests/test_minidsp.c`** — New test functions added to existing sections and registered in `main()`.
- No library API changes, no new files, no dependency changes.
