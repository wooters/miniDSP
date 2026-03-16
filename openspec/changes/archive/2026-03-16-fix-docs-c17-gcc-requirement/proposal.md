## Why

The codebase was downgraded from C23 to C17 (commit `9933eca`) for broader CI compatibility, but the README still references C23 and states that GCC 14 is required on Ubuntu. This misleads users into thinking they need a newer compiler than what's actually needed. C17 is supported by GCC 8+ and the default GCC 13 on Ubuntu 24.04, so the GCC 14 requirement no longer applies.

## What Changes

- Update the language badge from "C23" to "C17"
- Change `-std=c23` to `-std=c17` in the "Use in your project" compile examples
- Remove or rewrite the paragraph about needing GCC 14 on Ubuntu
- Update the container-test description that mentions "GCC 14"
- Simplify the Dockerfile to use the system default `gcc` instead of explicitly installing `gcc-14`

## Capabilities

### New Capabilities

(none — this is a docs/config fix only)

### Modified Capabilities

(none — no spec-level behavior changes)

## Impact

- `README.md` — multiple sections updated
- `Dockerfile` — simplified to drop `gcc-14` explicit install and `update-alternatives`
- `CLAUDE.md` — update the "C17 Notes" if it references C23 anywhere
- No library code or API changes
