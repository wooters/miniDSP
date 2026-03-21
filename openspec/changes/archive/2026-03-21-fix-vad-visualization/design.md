## Context

The VAD guide (`guides/vad.md`) uses the page anchor `{#vad}`, causing Doxygen to emit the rendered page as `vad.html`. The example program `examples/vad.c` also generates a Plotly visualization named `vad.html`, listed in `HTML_EXTRA_FILES`. Both files land in `docs/html/` — one overwrites the other, and the iframe ends up embedding the guide page recursively.

## Goals / Non-Goals

**Goals:**
- Eliminate the filename collision so the iframe shows the interactive visualization
- Minimal, localized change — rename the visualization file only

**Non-Goals:**
- Changing the Doxygen page anchor (`{#vad}`) — this would break any existing bookmarks or cross-references
- Refactoring the VAD example program beyond the output filename

## Decisions

**Rename the visualization file to `vad_plot.html`**

This follows the existing naming pattern used by other examples (e.g., the convention of `<name>_plot.html` for visualization assets distinct from guide pages). Renaming the visualization rather than the page anchor avoids breaking Doxygen cross-references (`\ref vad`, `\subpage vad`).

Alternative considered: Rename the page anchor to `{#vad-guide}`. Rejected because it would change the Doxygen-generated URL, potentially breaking existing links and requiring updates to `\subpage` references in `guides/tutorials.md`.

## Risks / Trade-offs

- [Low risk] Anyone with a bookmark to the old visualization URL will get a 404 → Acceptable since the visualization is embedded in the guide, not linked directly.
