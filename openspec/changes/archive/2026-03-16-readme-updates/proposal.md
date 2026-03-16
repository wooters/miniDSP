## Why

The README has several inaccuracies and omissions: it doesn't mention supported platforms, is missing the `audio_steg` tool, has sections in a suboptimal order (usage before build instructions), and incorrectly states the Apple `container` CLI is "built-in" on macOS 26+.

## What Changes

- Add a platform compatibility note near the top (Ubuntu, macOS; Windows untested — PRs welcome)
- Add `audio_steg` tool to the "Tools" section alongside `mel_viz`
- Move "Build and Test" section above "Use in your project"
- Fix the Apple `container` dependency row: it must be installed from GitHub, not built-in

## Capabilities

### New Capabilities

_(none — this is a documentation-only change)_

### Modified Capabilities

_(none — no spec-level behavior changes)_

## Impact

- `README.md` only — no code, API, or build system changes
