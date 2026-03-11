## Why

The mel_viz tool produces compelling real-time visualizations synced to audio, but there's no way to share them outside a browser. Users want to export the visualization as a video file (with audio) for sharing on social media, embedding in presentations, or archiving. This was an intended feature that was never implemented.

## What Changes

- Add an **Export** button to the mel_viz web UI (file mode only) that renders the visualization offline (faster than real-time) and produces a downloadable MP4 video file with audio.
- Use the **WebCodecs API** (`VideoEncoder` + `AudioEncoder`) with the **mp4-muxer** library to encode and mux frames offline — no real-time playback required.
- The renderer steps through precomputed mel frames programmatically, capturing each canvas frame and encoding it. Audio is decoded separately via `AudioContext.decodeAudioData()`.
- A progress bar shows export progress. The user can cancel mid-export.
- A small note near the Export button indicates Chrome/Edge-only support (WebCodecs requirement).

## Capabilities

### New Capabilities
- `video-export`: Offline rendering of the canvas visualization with audio into a downloadable MP4 video file using WebCodecs and mp4-muxer.

### Modified Capabilities
<!-- No existing spec-level behavior changes -->

## Impact

- **Files modified**: `tools/mel_viz/web/index.html` (export button, orchestration), `tools/mel_viz/web/style.css` (export UI styling). New `tools/mel_viz/web/exporter.js` module.
- **APIs**: No C library changes. Pure frontend/browser feature.
- **Dependencies**: [mp4-muxer](https://github.com/Vani-GitHub/mp4-muxer) (~30KB) — a lightweight JS library for muxing WebCodecs output into MP4 containers. Loaded via CDN or bundled.
- **Browser support**: Chrome 94+ and Edge 94+ (WebCodecs). Not supported in Firefox or Safari. Export button hidden in unsupported browsers with a note explaining the limitation.
