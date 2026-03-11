## Context

The mel_viz tool has a two-part architecture: a C backend (`mel_viz.c`) that precomputes mel energies from a WAV file into `data.js`, and a web frontend that renders a Canvas 2D radial visualization synced to an `<audio>` element. The visualization loop uses `requestAnimationFrame` and reads mel frame data from the precomputed array based on `audio.currentTime`.

There is currently no way to capture or export the visualization. Users must use external screen-recording tools, which produce inferior quality (frame drops, cursor artifacts, system audio mixing).

The key architectural insight: since file mode uses **precomputed mel data**, the renderer doesn't need real-time audio playback. It just needs a frame index. This enables offline (faster-than-real-time) export.

## Goals / Non-Goals

**Goals:**
- One-click export of the canvas visualization with synchronized audio as a downloadable MP4 video file.
- Offline rendering — export runs faster than real-time, no need to play through the entire audio.
- Every frame guaranteed (no drops regardless of system load).
- Progress indication and cancel support during export.

**Non-Goals:**
- Export from mic mode (no precomputed data to replay deterministically).
- Firefox/Safari support (WebCodecs not available — future enhancement via FFmpeg.wasm).
- Configurable resolution, bitrate, or codec selection.
- Server-side rendering or transcoding.

## Decisions

### 1. Encoding API: WebCodecs (VideoEncoder + AudioEncoder)

**Choice**: Use the WebCodecs API to encode canvas frames and audio samples offline, then mux with mp4-muxer.

**Rationale**: Enables faster-than-real-time export. The renderer steps through precomputed frames programmatically — no playback required. Every frame is guaranteed (no drops). Hardware-accelerated encoding where available.

**Alternatives considered**:
- *MediaRecorder + canvas.captureStream()*: Simpler, broader browser support, but requires real-time playback (30s song = 30s export). Unacceptable UX for longer audio.
- *FFmpeg.wasm*: Broad browser support, but ~25MB download. Overkill for v1 of a dev tool.

### 2. Output format: MP4 (H.264 video + AAC audio)

**Choice**: Export as MP4 with H.264 video and AAC audio via mp4-muxer.

**Rationale**: MP4 is universally playable (every OS, device, and social media platform). H.264 encoding is hardware-accelerated on most systems. mp4-muxer is a lightweight (~30KB) purpose-built library for muxing WebCodecs output.

**Alternative considered**:
- *WebM (VP8/VP9 + Opus)*: Less universally playable. webm-muxer is similar in size, but MP4 is the better choice for shareability.

### 3. Offline render loop

**Choice**: Decouple the renderer from `requestAnimationFrame`. In export mode, loop through frame indices 0..N synchronously (or in batched microtasks), calling `renderFrame(frameIndex / fps)` for each frame, then encoding the canvas contents via `VideoEncoder`.

**Rationale**: The existing `renderFrame()` already accepts a time parameter and indexes into the precomputed `MEL_VIZ_DATA.frames` array. The only change needed is feeding it simulated timestamps instead of `audio.currentTime`.

**Implementation detail**: The EMA smoothing state in the renderer must be reset before export starts so the exported video matches a fresh playback from the beginning.

### 4. Audio handling: Decode WAV separately

**Choice**: Fetch `audio.wav` via `fetch()`, decode with `AudioContext.decodeAudioData()` to get raw PCM (`AudioBuffer`), then encode PCM chunks with `AudioEncoder` (AAC codec).

**Rationale**: Decouples audio from the `<audio>` element entirely. The raw PCM can be fed to `AudioEncoder` in chunks without any playback.

### 5. Export module: Separate `exporter.js` file

**Choice**: Create a new `tools/mel_viz/web/exporter.js` module encapsulating all export logic.

**Rationale**: Keeps `index.html` focused on orchestration. The exporter manages its own state (idle/exporting/done) and exposes `startExport(canvas, renderFrameFn, melData, audioUrl, onProgress)` / `cancelExport()` methods.

### 6. Browser compatibility note

**Choice**: Display a small text note near the Export button: "Export requires Chrome or Edge".

**Rationale**: WebCodecs is Chrome/Edge-only. Rather than silently hiding the feature, a visible note sets expectations. The Export button itself is still hidden in unsupported browsers.

## Risks / Trade-offs

- **[Browser support]** WebCodecs is Chrome/Edge only (no Firefox, no Safari). → Mitigation: Feature-detect at runtime; hide Export button and show note in unsupported browsers. Acceptable for a dev tool. FFmpeg.wasm can be added as a fallback later.
- **[mp4-muxer dependency]** Adds a third-party JS dependency (~30KB). → Mitigation: Well-maintained library, small footprint, loaded via CDN with integrity hash. No build tooling required.
- **[EMA smoothing state]** The renderer's temporal smoothing accumulates state. Export must reset this state to match a fresh playback. → Mitigation: Add a `resetSmoothing()` function to the renderer, called before export begins.
- **[Canvas resolution]** Export captures at the current canvas resolution. If the canvas is small (e.g., on a narrow viewport), the video will be low-res. → Mitigation: Acceptable for v1. Could add a resolution override in the future.
- **[Memory usage]** For long audio files, holding all encoded chunks in memory before muxing could use significant RAM. → Mitigation: mp4-muxer supports streaming mode (writes to an ArrayBuffer incrementally). Use this for large files.
