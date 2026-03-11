## Why

The `examples/` directory mixes simple library demonstrations (sine waves, impulse responses) with substantial standalone programs that have multi-mode CLIs, embedded assets, and independent use cases. This blurs the distinction between "learn the API" examples and "do something useful" tools. `audio_steg.c` (641 lines, 6 operational modes, embedded PNG asset) is clearly a tool — not a teaching example. Moving it to `tools/` clarifies both directories and follows the pattern established by `tools/mel_viz/`.

## What Changes

- Move `examples/audio_steg.c` and its embedded asset (`examples/space_invader.png`) to `tools/audio_steg/` with its own Makefile and README
- Remove `audio_steg` from `examples/Makefile` (EXAMPLES variable, plot target)
- Remove `audio_steg` entries from `.gitignore` and `.dockerignore`
- Add `audio_steg` build to root Makefile `tools` target
- Update Doxygen configuration if `audio_steg.c` is referenced in `EXAMPLE_PATH` or guides
- Update `README.md` to reflect the move

## Capabilities

### New Capabilities
- `tool-audio-steg`: Relocate audio_steg from examples to tools with its own build system, README, and directory structure following the mel_viz pattern.

### Modified Capabilities

_(none — no spec-level behavior changes to existing capabilities)_

## Impact

- **Build system**: `examples/Makefile` loses one target; `tools/audio_steg/Makefile` is added; root Makefile `tools` target updated
- **File layout**: `examples/audio_steg.c` and `examples/space_invader.png` move to `tools/audio_steg/`
- **Documentation**: README, Doxygen guides, and any snippet references pointing at `examples/audio_steg.c` need path updates
- **Git tracking**: `.gitignore` and `.dockerignore` entries change from `examples/audio_steg*` to `tools/audio_steg/` patterns
- **No API changes**: The library itself is untouched
