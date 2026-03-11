## 1. Dependencies & Setup

- [x] 1.1 Add mp4-muxer to the project (CDN `<script>` tag in `index.html` or vendored copy in `tools/mel_viz/web/`) — Used Mediabunny (same author, higher-level API) via esm.sh CDN import in exporter.js
- [x] 1.2 Create `tools/mel_viz/web/exporter.js` with WebCodecs feature detection (`typeof VideoEncoder !== 'undefined'` and `typeof AudioEncoder !== 'undefined'`)

## 2. Offline Render Loop

- [x] 2.1 Add a `resetSmoothing()` function to `renderer.js` that clears the EMA smoothing state — Existing `resetState()` already does this; reused as-is
- [x] 2.2 Implement the frame-stepping loop in `exporter.js`: iterate frame indices 0..N, call `renderFrame(frameIndex / fps)` for each, yield to the event loop periodically (e.g., every 10 frames via `setTimeout(0)`) to keep the UI responsive
- [x] 2.3 After each `renderFrame()`, create a `VideoFrame` from the canvas and pass it to `VideoEncoder.encode()` — Handled by Mediabunny's CanvasSource.add()

## 3. Audio Encoding

- [x] 3.1 Fetch `audio.wav` and decode with `AudioContext.decodeAudioData()` to get raw PCM (`AudioBuffer`)
- [x] 3.2 Convert `AudioBuffer` channel data into `AudioData` objects and encode with `AudioEncoder` (AAC codec) — Handled by Mediabunny's AudioBufferSource.add()

## 4. Muxing & Download

- [x] 4.1 Initialize `mp4-muxer` `Muxer` with H.264 video track and AAC audio track — Using Mediabunny Output with Mp4OutputFormat
- [x] 4.2 Feed encoded video and audio chunks to the muxer as they arrive from the encoders
- [x] 4.3 On completion, finalize the muxer, create a Blob from the output, and trigger download of `mel-viz-export.mp4`

## 5. UI Integration

- [x] 5.1 Add `<script>` tags for mp4-muxer and `exporter.js` in `index.html` — ES module import in exporter.js handles the CDN dependency
- [x] 5.2 Add Export button to the controls area, conditionally shown only in file mode when WebCodecs is available
- [x] 5.3 Add browser compatibility note near the Export button — Updated to "Not supported in Firefox" (Safari works too)
- [x] 5.4 Toggle button label between "Export" and "Cancel" based on export state
- [x] 5.5 Add progress bar showing `framesProcessed / totalFrames` during export
- [x] 5.6 Wire up Cancel button to abort encoders, discard output, and reset UI state

## 6. Styling

- [x] 6.1 Add CSS styles for the Export/Cancel button in `style.css`
- [x] 6.2 Add CSS for the progress bar and browser compatibility note

## 7. Testing & Verification

- [x] 7.1 Manual test: open mel_viz with a WAV file in Chrome, click Export, verify MP4 downloads after export completes
- [x] 7.2 Manual test: click Cancel mid-export, verify no download and UI resets
- [x] 7.3 Manual test: verify Export button does not appear in mic mode
- [x] 7.4 Manual test: verify Export button does not appear in Firefox/Safari — Safari supports WebCodecs; button appears and export works. Firefox not tested (unavailable).
- [x] 7.5 Manual test: verify exported MP4 plays correctly in VLC, QuickTime, and browser
- [x] 7.6 Manual test: verify progress bar advances during export
