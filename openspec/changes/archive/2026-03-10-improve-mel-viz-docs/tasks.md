## 1. Media Assets

- [x] 1.1 Create `docs/assets/` directory
- [x] 1.2 Compress video with ffmpeg: scale to 480x480, CRF 33, baseline profile, AAC 64k audio, faststart — output to `docs/assets/mel-viz-demo.mp4` (3.5 MB)
- [x] 1.3 Take screenshot of mel_viz web UI via Puppeteer and save to `docs/assets/mel-viz-screenshot.png`
- [x] 1.4 Add both files to `Doxyfile` `HTML_EXTRA_FILES`

## 2. Guide Page Rewrite

- [x] 2.1 Rewrite `guides/mel-viz.md` with new structure: title, video embed (`\htmlonly` + `<video controls loop playsinline>`), screenshot (`\image html`), tutorial, mic mode, architecture
- [x] 2.2 Write "Visualize your own audio" tutorial section with numbered steps: build, run on WAV, serve output, customize controls, export video
- [x] 2.3 Move architecture diagram and control reference to end of page

## 3. Verification

- [x] 3.1 Run `doxygen` and confirm video and screenshot appear in `docs/html/`
- [x] 3.2 Open generated guide page in browser and verify video plays and screenshot renders
