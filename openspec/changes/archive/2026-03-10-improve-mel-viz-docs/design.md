## Context

The mel_viz guide (`guides/mel-viz.md`) currently reads like a reference card: build command, file mode snippet, mic mode snippet, architecture diagram, one-liner about controls. There are no visuals showing what the tool actually produces. The user has an exported MP4 demo video (13 MB, 30s, 1080x1080 H.264) and wants a screenshot and a step-by-step tutorial added.

GitHub has a 100 MB hard file limit and warns above 50 MB, but 13 MB is still heavy for a docs asset in a small C library repo. Compression is warranted.

## Goals / Non-Goals

**Goals:**
- Show what mel_viz looks like before the reader builds anything (video + screenshot)
- Provide a complete, numbered walkthrough for going from "I have a WAV file" to "I see a visualization in my browser"
- Keep media assets small enough to commit without LFS

**Non-Goals:**
- Rewriting the tool's standalone README (`tools/mel_viz/README.md`) — keep it as the CLI/architecture reference
- Adding video/screenshot to the README (Doxygen guide only)
- Hosting media externally (keep everything in-repo for simplicity)

## Decisions

### 1. Compress video with ffmpeg before committing

**Decision**: Re-encode the MP4 to ~3.5 MB using ffmpeg: scale to 480x480, CRF 33, H.264 baseline profile, AAC audio at 64 kbps.

**Rationale**: 13 MB is fine for GitHub but disproportionate for a docs asset in a ~2 MB repo. 480px is plenty for an inline demo video. CRF 33 is visually acceptable for a colorful animation (no text to read). Baseline profile ensures maximum browser compatibility.

**Command**:
```sh
ffmpeg -i ~/Downloads/mel-viz-export.mp4 \
  -vf scale=480:480 -c:v libx264 -crf 33 -profile:v baseline \
  -c:a aac -b:a 64k -movflags +faststart docs/assets/mel-viz-demo.mp4
```

Audio is kept (AAC 64 kbps) so the viewer can hear the music driving the visualization. `-movflags +faststart` enables progressive playback.

**Alternatives considered**:
- Git LFS: adds infrastructure complexity for one file
- External hosting (YouTube/S3): breaks self-contained repo, links rot
- GIF: much larger file for same quality, no playback controls

### 2. Screenshot via Puppeteer during implementation

**Decision**: Take a browser screenshot of the running mel_viz visualization using Puppeteer and save as PNG. Crop/resize to ~800px wide.

**Rationale**: A real screenshot of the actual tool is more authentic than a mockup. Puppeteer is already available as an MCP tool.

**Destination**: `docs/assets/mel-viz-screenshot.png`

### 3. Media stored in `docs/assets/`

**Decision**: Create `docs/assets/` directory for documentation media files. Add both files to `Doxyfile` `HTML_EXTRA_FILES` so Doxygen copies them to the output directory.

**Rationale**: Keeps media separate from source code. `docs/` is already used for Doxygen output (`docs/html/`). The `assets/` subdirectory makes it clear these are source assets, not generated output.

### 4. Guide page structure: visuals first, tutorial second, reference last

**Decision**: Restructure `guides/mel-viz.md` as:
1. One-line description
2. Video embed (HTML5 `<video>`, autoplay muted loop)
3. Screenshot (for non-video-capable contexts)
4. "Visualize your own audio" tutorial (numbered steps)
5. Live mic mode (kept brief)
6. Architecture diagram (moved to end)

**Rationale**: Lead with the "wow" — show what it does before explaining how. The tutorial is the primary call-to-action. Architecture is reference material for those who want to understand internals.

### 5. Video embedded via `\htmlonly` + HTML5 `<video>`

**Decision**: Use Doxygen's `\htmlonly`/`\endhtmlonly` with an HTML5 `<video>` tag. Attributes: `controls loop playsinline width="480"`. No autoplay — the video has audio, so the user should click play.

**Rationale**: Doxygen has no built-in video command. `\htmlonly` is the established pattern in this project (see MEMORY.md notes on HTML5 media). `autoplay muted` is required by browsers for autoplay. `controls` lets users pause/scrub. `loop` keeps the demo running.

### 6. Screenshot embedded via standard Doxygen `\image`

**Decision**: Use `\image html mel-viz-screenshot.png "mel_viz visualization" width=600px` for the screenshot.

**Rationale**: Unlike `<video>`, Doxygen natively supports `\image` for static images. Simpler than `\htmlonly`.

## Risks / Trade-offs

- **[Video codec compatibility]** → H.264 baseline profile works in all modern browsers. Safari and Firefox both support it. Low risk.
- **[Repo size growth]** → ~2-3 MB video + ~200 KB screenshot adds ~3 MB to repo. Acceptable for a one-time addition. If more media assets accumulate in the future, consider LFS then.
- **[Screenshot staleness]** → If the UI changes, the screenshot becomes outdated. Mitigation: screenshot is easy to retake; note the Puppeteer approach so future maintainers can redo it.
- **[Audio in video]** → Audio is included at 64 kbps AAC so the viewer can hear the music driving the visualization. Adds ~240 KB to the file, well worth it for the demo experience.
