## 1. Project scaffolding

- [x] 1.1 Create `tools/mel_viz/` directory structure (`mel_viz.c`, `Makefile`, `web/`)
- [x] 1.2 Create `tools/mel_viz/Makefile` (build mel_viz binary, link libminidsp + libsndfile + fftw3)
- [x] 1.3 Add `tools` target to root Makefile that delegates to `tools/mel_viz/Makefile`
- [x] 1.4 Update `.gitignore` and `.dockerignore` for `tools/mel_viz/mel_viz` binary and output artifacts

## 2. C program — WAV analysis and output assembly

- [x] 2.1 Implement CLI argument parsing (input WAV, `-o`, `--mels`, `--fft-size`, `--fps`, `--min-freq`, `--max-freq`, `--width`, `--height`)
- [x] 2.2 Read WAV file via libsndfile (mono downmix if stereo, read sample rate)
- [x] 2.3 Compute mel energies per video frame: slide window by `sample_rate / fps` samples, call `MD_mel_energies()` per frame
- [x] 2.4 Compute bass envelope: sum lowest mel bands per frame, apply one-pole envelope follower
- [x] 2.5 Write `data.js` — emit `MEL_VIZ_DATA` object with frames array, bass envelope, metadata
- [x] 2.6 Assemble output directory: create dir, copy `web/*` assets, copy input WAV as `audio.wav`, write `data.js`
- [x] 2.7 Print summary (num frames, duration, output path) and usage hint (`open <dir>/index.html`)

## 3. Web renderer — HTML scaffold and styles

- [x] 3.1 Create `web/index.html` — canvas element, audio element, control panel, script module imports
- [x] 3.2 Create `web/style.css` — dark layout, canvas centering, control panel styling, responsive sizing

## 4. Web renderer — Canvas 2D radial engine

- [x] 4.1 Create `web/renderer.js` — main render loop with `requestAnimationFrame`
- [x] 4.2 Implement background layer: radial gradient with hue/brightness driven by total energy
- [x] 4.3 Implement ring drawing: concentric rings with radius = base + energy * scale, ring thickness proportional to energy
- [x] 4.4 Implement glow effect: `ctx.shadowBlur` and `ctx.shadowColor` per ring, intensity driven by energy
- [x] 4.5 Implement ring wobble: sinusoidal radius deformation (`radius + sin(angle * n) * wobble * energy`)
- [x] 4.6 Implement center bloom: bright circle at center, radius and opacity driven by bass energy
- [x] 4.7 Implement bass pulse: global `ctx.scale()` driven by bass envelope (attack/decay smoothed)
- [x] 4.8 Implement temporal smoothing: EMA filter on mel energies, controlled by smoothing knob
- [x] 4.9 Implement mel band grouping: average adjacent mel bands into visual ring groups (8/12/24 presets)

## 5. Web renderer — Data providers

- [x] 5.1 Create `web/audio-provider.js` — `MelFrameProvider` interface with `getFrame(time)` method
- [x] 5.2 Implement file mode: index into precomputed `MEL_VIZ_DATA.frames` array by `audio.currentTime * fps`
- [x] 5.3 Implement mic mode: `getUserMedia()` + `AnalyserNode` + `getFloatFrequencyData()` + JS mel weighting
- [x] 5.4 Auto-detect mode: if `MEL_VIZ_DATA` exists → file mode with mic toggle; otherwise → mic-only mode

## 6. Web renderer — Palettes and controls

- [x] 6.1 Create `web/palettes.js` — define plasma, ocean, fire, neon palettes (bgHue, ring colors, glow color)
- [x] 6.2 Create `web/controls.js` — slider UI for smoothing, bass sensitivity, wobble, glow intensity
- [x] 6.3 Wire palette dropdown and ring-group dropdown to renderer
- [x] 6.4 Wire sliders to renderer parameters with real-time updates (no recomputation)

## 7. Integration and polish

- [x] 7.1 End-to-end test: run `mel_viz` on a test WAV, open output, verify audio-synced playback
- [x] 7.2 Test live mic mode: open `web/index.html` directly (via local server), verify mic visualization
- [x] 7.3 Test all four palettes and all knob ranges
- [x] 7.4 Verify seeking/scrubbing audio timeline syncs visual correctly
- [x] 7.5 Add `tools/mel_viz/README.md` with build/usage instructions and screenshots

## 8. Build system integration

- [x] 8.1 Verify `make tools` builds `mel_viz` from a clean state
- [x] 8.2 Verify `make clean` (or `make -C tools/mel_viz clean`) removes build artifacts
- [x] 8.3 Test on macOS (Apple Silicon + Homebrew paths) and Linux
