## Why

The VAD guide page (`guides/vad.md`) uses the Doxygen page anchor `{#vad}`, which causes Doxygen to emit the rendered guide as `vad.html`. The `HTML_EXTRA_FILES` entry also copies the example-generated `examples/vad.html` (the interactive Plotly visualization) into the same output directory. One file overwrites the other, so the iframe `<iframe src="vad.html">` ends up embedding the guide page itself instead of the visualization — a recursive embed.

## What Changes

- Rename the visualization output file from `vad.html` to `vad_plot.html` to eliminate the filename collision with Doxygen's generated `vad.html` page.
- Update all references: the C example program, Doxyfile `HTML_EXTRA_FILES`, the guide iframe `src`, `.gitignore`, and `.dockerignore`.

## Capabilities

### New Capabilities

(none)

### Modified Capabilities

(none — this is a bug fix with no spec-level behavior changes)

## Impact

- `examples/vad.c` — change output filename from `vad.html` to `vad_plot.html`
- `guides/vad.md` — update iframe `src` attribute
- `Doxyfile` — update `HTML_EXTRA_FILES` entry
- `.gitignore` / `.dockerignore` — update ignore patterns
