## Why

The mel_viz tool is the first entry in the new `tools/` directory but has no presence in the project's documentation. The Doxygen-generated docs site has tutorials for every example, and the README lists all capabilities, but neither mentions mel_viz. Additionally, there's no "Tools" section in either place — since this is the first tool (as opposed to an example), we need to establish the convention.

## What Changes

- Add a Doxygen guide page for mel_viz (`guides/mel-viz.md`) explaining the tool, its architecture, and usage
- Add a "Tools" section to `guides/tutorials.md` (separate from the existing tutorial list) with a `\subpage` link
- Add a "Tools" section to `README.md` listing mel_viz with a brief description and usage snippet
- Update `Doxyfile` if needed (INPUT paths, EXAMPLE_PATH, etc.)

## Capabilities

### New Capabilities
- `mel-viz-docs`: Doxygen guide page, README section, and tutorials index entry for the mel_viz audio visualizer tool

### Modified Capabilities
<!-- None -->

## Impact

- **Documentation**: New guide page, updated tutorials index, updated README
- **Doxyfile**: May need updates if guide references assets or snippets from `tools/`
- **No code changes**: This is purely documentation
