## Context

A full audit of 206+ division operations across 11 source files found two unguarded division-by-zero bugs reachable with valid API inputs:

1. **`MD_adjust_dblevel`** (`src/minidsp_core.c:177`): `gain = sqrt(desired_power * N / input_energy)` where `input_energy == 0.0` for silent input.
2. **`BiQuad_new`** (`src/biquad.c:113`): `alpha = sn * sinh(M_LN2/2 * bandwidth * omega / sn)` where `sn = sin(omega) == 0.0` when `freq == 0` or `freq == srate/2`.

All other divisions are safe (guarded by asserts, conditionals, epsilon offsets, or mathematical guarantees).

## Goals / Non-Goals

**Goals:**
- Guard both functions against the identified zero-denominator conditions.
- Produce sensible output (not NaN/Inf/garbage) for the triggering inputs.
- Add tests that verify the guarded behavior.

**Non-Goals:**
- Not adding guards to already-safe divisions (would be noise).
- Not changing any function signatures or return types.
- Not auditing third-party code (biquad is project-owned).

## Decisions

**`MD_adjust_dblevel` with zero energy: copy input unchanged.**
When `input_energy == 0.0`, the function copies `in` to `out` and returns. Rationale: you cannot amplify silence to any target dB — the output must remain silence. This is the only mathematically correct result.

*Alternative considered:* Return an error code. Rejected — the function returns `void`, and changing to `int` would be a breaking API change for a rare edge case. The "copy silence" behavior is correct, not an error.

**`BiQuad_new` with freq=0 or freq=Nyquist: return NULL.**
When `sn == 0.0` (i.e., `freq == 0` or `freq == srate/2`), the function returns `NULL` (after freeing the allocated struct). This is consistent with the existing error path for invalid filter types.

*Alternative considered:* Clamp frequency to `[1, srate/2 - 1]`. Rejected — silently modifying the user's requested frequency would violate least-surprise. `NULL` signals the caller that the parameters are invalid, matching the existing convention.

**Guard placement: early return before the division, not epsilon offset.**
An `if (denominator == 0.0)` check before the division is clearer than adding `DBL_MIN` to the denominator. Epsilon offsets produce technically finite but astronomically large/wrong results. An early return produces the correct behavior.

## Risks / Trade-offs

**Risk: Existing callers of `MD_adjust_dblevel` pass silent signals.** → The new behavior (output = silence) is strictly better than the old behavior (output = NaN). No caller could have been relying on NaN output.

**Risk: Existing callers of `BiQuad_new` pass freq=0.** → Previously this produced garbage coefficients that would corrupt audio. Returning NULL is strictly better. However, callers that don't check the return value will segfault. The existing `test_biquad_invalid_type` test already demonstrates the NULL-return convention, so this is documented behavior.
