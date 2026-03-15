## Why

The documentation currently assumes a student audience ("so students can start with…", "for students"). miniDSP is useful to anyone doing audio DSP — hobbyists, engineers, researchers — not just students. Replacing "student" with "user" makes the docs more inclusive without changing their substance.

## What Changes

- Replace "students" → "users" and "for students" → "for users" in all documentation and source comments.
- Two occurrences total:
  - `guides/fir-convolution.md` line 6: "so students can start with…"
  - `src/biquad.c` line 17: "How a biquad filter works (for students):"

## Capabilities

### New Capabilities

(none)

### Modified Capabilities

(none — this is a wording-only change with no behavioral or requirement impact)

## Impact

- Two files changed: `guides/fir-convolution.md`, `src/biquad.c`
- No API, build, or test impact
- Regenerated Doxygen output will reflect the updated wording
