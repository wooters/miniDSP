## Why

`MD_lowpass_brickwall()` was added to the library but has no tutorial-style guide page. It is documented only in the header doc-comment and mentioned in passing within the audio steganography pipeline diagram. Users browsing the tutorials page have no way to discover or learn about this function.

## What Changes

- **New guide page** for the brickwall lowpass filter, added to the FIR/filtering section of the tutorials index. Covers the theory (FFT-based frequency-domain zeroing), the formula with a "Reading the formula in C" section, API usage, Gibbs ringing tradeoffs, and a working example program with before/after spectrum visualization.
- **Tutorials index updated** — add a `\subpage` entry for the new brickwall guide in `guides/tutorials.md`, positioned near the existing FIR convolution and resampling entries.

## Capabilities

### New Capabilities
- `brickwall-docs`: Tutorial guide page for `MD_lowpass_brickwall()` — theory, formula, C walkthrough, API, example program, and integration into the tutorials index.

### Modified Capabilities
_(none — no existing spec requirements change)_

## Impact

- New file: `guides/brickwall-lowpass.md`
- Modified file: `guides/tutorials.md` (add `\subpage` entry)
- New example program: `examples/brickwall.c` with snippet markers for guide embedding
- Modified file: `examples/Makefile` (add to `EXAMPLES`, `plot` target)
- Modified files: `.gitignore`, `.dockerignore` (add example artifact lines)
- No API or library code changes.
