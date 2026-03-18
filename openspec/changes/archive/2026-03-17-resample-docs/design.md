## Context

The `mel_viz` and `audio_steg` tools each have a Doxygen guide page under `guides/` and a `\subpage` entry in `guides/tutorials.md`. The new `resample` tool has no guide page yet.

## Goals / Non-Goals

**Goals:**
- Add a guide page following the same structure as `mel-viz.md` (the simpler of the two tool guides)
- Document usage, options, example workflows, and a brief explanation of how the resampler works
- Make the page discoverable from the Tutorials navigation

**Non-Goals:**
- Deep mathematical treatment of sinc interpolation (already covered in `guides/resampling.md`)
- Generated audio assets or interactive plots (this is a CLI tool, not a signal generator)
- Adding the tool to CI docs builds (it's already built by `make tools`)

## Decisions

**1. Guide structure**: Follow the `mel-viz.md` pattern — short, practical, focused on running the tool. Link to `guides/resampling.md` for the underlying API and math.

*Rationale*: The mel-viz guide is a clean tool-documentation template. The audio-steganography guide mixes API docs with tool docs, which is appropriate for that tool but overkill here since the resampling API already has its own guide.

**2. No generated assets**: The resample tool produces WAV files, not visualizations. No audio samples or plots to embed.

## Risks / Trade-offs

None — documentation-only change with no code impact.
