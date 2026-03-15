## Why

A full audit of all 206+ division operations across 11 source files found two unguarded division-by-zero bugs that can be triggered with valid (non-assert-guarded) inputs. `MD_adjust_dblevel` produces NaN/Inf when given a silent signal, and `BiQuad_new` divides by `sin(omega)` which is zero when the filter frequency is 0 Hz or exactly Nyquist. Both are reachable through normal API usage without violating documented preconditions.

## What Changes

- **Fix `MD_adjust_dblevel` (HIGH):** Guard against zero input energy. When the input signal is all zeros (or near-zero), the function currently computes `gain = sqrt(desired_power * N / 0.0)`, producing NaN. The fix will detect zero energy and copy the input to output unchanged (you can't amplify silence to a target dB level).
- **Fix `BiQuad_new` (MEDIUM):** Guard against `sin(omega) == 0` in the biquad coefficient calculation. This occurs when `freq == 0` or `freq == srate/2`. The fix will return `NULL` for these degenerate filter frequencies, consistent with the existing invalid-type error path.
- **Add tests** for both fixes to `tests/test_minidsp.c`.

## Capabilities

### New Capabilities
- `division-by-zero-guards`: Guards against division-by-zero in library functions that accept valid inputs which can produce zero denominators

### Modified Capabilities

## Impact

- **`src/minidsp_core.c`** — `MD_adjust_dblevel()`: add zero-energy guard
- **`src/biquad.c`** — `BiQuad_new()`: add zero-sin(omega) guard
- **`tests/test_minidsp.c`** — new tests for both edge cases
- **No API signature changes** — return types and parameters are unchanged
- **Behavioral change:** `MD_adjust_dblevel` with silent input will now produce silent output (instead of NaN). `BiQuad_new` with freq=0 or freq=Nyquist will now return NULL (instead of producing garbage coefficients).
