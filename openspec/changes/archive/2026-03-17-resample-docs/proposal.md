## Why

The `mel_viz` and `audio_steg` tools both have Doxygen guide pages, but the new `resample` tool does not. Adding a guide page keeps the documentation consistent and makes the tool discoverable from the Tutorials/Tools navigation.

## What Changes

- Add a new Doxygen guide page `guides/resample-tool.md` documenting the resample CLI tool (usage, options, examples, how it works).
- Link it from `guides/tutorials.md` under the Tools section.
- Add `tools/resample` to the Doxyfile `EXAMPLE_PATH` so `\snippet` directives can reference `resample.c`.

## Capabilities

### New Capabilities
- `resample-tool-guide`: Doxygen guide page for the resample CLI tool — usage, options (`-z`, `-b`), example workflows, and how the polyphase sinc resampler works.

### Modified Capabilities
(none)

## Impact

- **New files**: `guides/resample-tool.md`
- **Modified files**: `guides/tutorials.md` (add `\subpage`), `Doxyfile` (`EXAMPLE_PATH`)
- **No code changes** — documentation only
