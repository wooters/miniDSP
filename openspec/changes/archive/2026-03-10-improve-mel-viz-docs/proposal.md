## Why

The mel_viz guide page is a dry reference listing — it tells you what mel_viz is and shows CLI flags, but doesn't show what the visualization actually looks like or walk someone through creating their own. A video demo, a screenshot, and a step-by-step tutorial would make the docs engaging and actionable.

## What Changes

- **Add video demo**: Embed the exported MP4 (30s, 1080x1080) showing the visualizer in action. The raw file is 13 MB — compress it to ~2-3 MB via ffmpeg (lower bitrate, smaller resolution) and commit directly (no LFS needed at that size). Host as a Doxygen `HTML_EXTRA_FILES` asset with an HTML5 `<video>` tag in the guide.
- **Add screenshot**: Capture a static screenshot of the mel_viz web UI and embed it in the guide. Gives immediate visual context even if the video doesn't autoplay.
- **Add step-by-step tutorial**: Write a "Visualize your own audio" section with numbered steps from WAV file to browser, including build, run, serve, and customize. Replace the current terse code blocks with a guided walkthrough.
- **Restructure the guide page**: Reorganize `guides/mel-viz.md` so it leads with visuals (video + screenshot), then the tutorial, then reference material (architecture, controls, CLI flags).

## Capabilities

### New Capabilities
- `mel-viz-media-assets`: Managing video and screenshot assets for mel_viz documentation (compression, Doxyfile integration, embedding)
- `mel-viz-tutorial`: Step-by-step tutorial for creating a visualization from a user's own audio file

### Modified Capabilities

_(none — no existing spec requirements are changing)_

## Impact

- `guides/mel-viz.md` — major rewrite
- `Doxyfile` — add compressed video and screenshot to `HTML_EXTRA_FILES`
- New media assets committed to repo (compressed video ~2-3 MB, screenshot ~200 KB)
- `tools/mel_viz/README.md` — may add cross-reference to the improved guide
