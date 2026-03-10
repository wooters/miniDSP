## Context

miniDSP is a C library for audio DSP. It already has `MD_mel_energies()` for per-frame mel-band analysis, `MD_stft()` for sliding-window FFT, and libsndfile-based audio I/O in examples. The visualizer combines these into a tool that produces browser-rendered audio visualizations.

## Goals / Non-Goals

**Goals:**
- Radial concentric-ring visualization driven by mel-band energies
- Reactive background that responds to overall energy
- Bass pulse effect (global zoom on kick/bass hits)
- Configurable color palettes (plasma, ocean, fire, neon)
- Interactive JS knobs for visual parameters (no recomputation needed)
- Live mic mode that works without the C program
- Configurable resolution and aspect ratio
- Folder-based output (HTML + JS + CSS + data + audio reference)

**Non-Goals:**
- WebGL rendering (Canvas 2D first; WebGL is a future upgrade path)
- npm/webpack/vite build tooling for the JS side (raw ES modules only)
- Streaming/real-time C-to-browser pipeline (the C program is batch; live mode is pure JS)
- Video file export from C (browser's MediaRecorder API handles this if needed)

## Decisions

### 1. Architecture: C computes, browser renders

**Decision**: The C program's sole job is DSP — read WAV, compute mel energies per video frame, write a `data.js` file, and copy web assets into an output folder. All rendering happens in the browser via Canvas 2D.

**Rationale**: Avoids C-side image rendering dependencies (stb_image, Cairo). Leverages the browser for animation, audio sync (`<audio>` + `currentTime`), interactivity, and cross-platform display. Already consistent with the project's existing HTML visualization pattern.

**Alternative considered**: C generates PNG frames + ffmpeg stitches to video. Rejected because it adds external dependencies, loses interactivity, and requires frame-accurate sync logic.

### 2. Output format: folder, not single file

**Decision**: Output is a directory containing `index.html`, JS modules, CSS, `data.js` (generated), and `audio.wav` (copied/symlinked).

**Rationale**: Keeps JS/CSS maintainable as separate files. The folder can be zipped for sharing. A `--embed` flag could inline everything into a single HTML in the future, but folder-first is cleaner for development and for the file sizes involved (~650KB data + multi-MB audio).

### 3. Mel band grouping

**Decision**: Compute 24 mel bands in C, group into 8 visual rings in JS. Grouping map is written into `data.js` and can be reconfigured via JS knobs.

**Mapping:**
| Ring | Mel Bands | Perceptual Region | Visual Role |
|------|-----------|-------------------|-------------|
| 0 | 0-2 | Sub-bass (20-100 Hz) | Core pulse |
| 1 | 3-5 | Bass (100-300 Hz) | Inner glow |
| 2 | 6-8 | Low-mid (300-800 Hz) | Body |
| 3 | 9-11 | Mid (800-2 kHz) | Warmth |
| 4 | 12-14 | Upper-mid (2-4 kHz) | Presence |
| 5 | 15-17 | Presence (4-7 kHz) | Brightness |
| 6 | 18-20 | Brilliance (7-12 kHz) | Shimmer |
| 7 | 21-23 | Air (12-16 kHz) | Outer halo |

**Rationale**: 24 raw rings is visually cluttered. 8 groups map well to perceptual regions and keep the visual clean. The grouping is done JS-side so presets (8, 12, 24 rings) can be offered as a knob.

### 4. Bass pulse effect

**Decision**: Sum mel bands 0-5 (rings 0-1) for a "bass energy" signal. Apply a one-pole envelope follower (`env = max(bass, env * decay)`) for smooth attack/fast-release. Use the envelope to scale the entire canvas via `ctx.scale()` — everything zooms toward the viewer on each kick.

**Rationale**: Simple, visually dramatic, and cheap to compute. The envelope follower prevents jitter while keeping the pulse snappy.

### 5. Live mic mode via Web Audio API

**Decision**: `web/index.html` is independently usable. When opened without `data.js`, it offers mic input via `getUserMedia()` + `AnalyserNode`. The JS applies mel weighting to FFT bins directly. The Canvas renderer doesn't know or care whether data comes from precomputed arrays or live audio.

**Rationale**: Makes the visualizer instantly demo-able without compiling C. The renderer uses a `MelFrameProvider` interface — file mode indexes into the precomputed array by audio `currentTime`; mic mode returns live-analyzed bins. Same renderer, two data sources.

### 6. Color palettes

**Decision**: Four built-in palettes, selectable via dropdown in the UI. Each palette defines: background hue, per-ring colors (array of 8), and glow color.

| Palette | Background | Bass → Air colors | Vibe |
|---------|-----------|-------------------|------|
| Plasma | Deep purple | Hot pink → orange → yellow → white | Lava lamp |
| Ocean | Deep blue | Navy → teal → cyan → ice white | Underwater |
| Fire | Dark red | Red → orange → yellow → white | Flames |
| Neon | Dark blue | Electric pink → cyan → purple → white | Synthwave |

**Rationale**: Palette as a knob was explicitly requested. Four gives good variety without being overwhelming. Adding more later is trivial (just add an object to `palettes.js`).

### 7. JS knobs (visual parameters)

**Decision**: Sliders/controls in the HTML UI for real-time visual tuning. No recomputation needed — these only affect the renderer.

| Knob | Type | Range | What it controls |
|------|------|-------|-----------------|
| Palette | Dropdown | plasma/ocean/fire/neon | Color scheme |
| Smoothing | Slider | 0.0-0.95 | Temporal EMA factor |
| Bass sensitivity | Slider | 0.0-1.0 | Pulse zoom depth |
| Ring wobble | Slider | 0.0-1.0 | Organic deformation |
| Glow intensity | Slider | 0.0-2.0 | Shadow blur multiplier |
| Ring groups | Dropdown | 8/12/24 | Number of visual rings |

### 8. C-side command-line interface

**Decision**: The C program accepts:
```
mel_viz <input.wav> [options]
  -o <dir>         Output directory (default: mel_viz_out/)
  --mels <n>       Number of mel bands (default: 24)
  --fft-size <n>   FFT window size (default: 2048)
  --fps <n>        Frames per second (default: 30)
  --min-freq <hz>  Low frequency bound (default: 40)
  --max-freq <hz>  High frequency bound (default: 16000)
  --width <px>     Canvas width (default: 1080)
  --height <px>    Canvas height (default: 1080)
```

### 9. Directory structure

```
tools/
└── mel_viz/
    ├── mel_viz.c              C program: WAV → data.js + output folder
    ├── Makefile               Builds mel_viz, copies web assets
    └── web/                   Browser renderer (pure JS, no build tools)
        ├── index.html         Main page: canvas + audio + controls
        ├── renderer.js        Canvas 2D radial drawing engine
        ├── audio-provider.js  File mode + mic mode data sources
        ├── palettes.js        Color palette definitions
        ├── controls.js        Slider/knob UI and state management
        └── style.css          Layout and control styling
```

### 10. `data.js` format

```javascript
const MEL_VIZ_DATA = {
  sampleRate: 44100,
  fps: 30,
  numFrames: 5400,
  numBands: 24,
  numGroups: 8,
  groupMap: [3, 3, 3, 3, 3, 3, 3, 3],
  width: 1080,
  height: 1080,
  // Flat array: frames[f * numBands + b] = energy
  frames: [0.001, 0.003, ...],
  // Per-frame bass energy (sum of bands 0-5, envelope-followed)
  bassEnvelope: [0.0, 0.001, ...],
  audioFile: "audio.wav"
};
```

### 11. Rendering layers (back to front)

1. **Background gradient** — radial gradient, hue from palette, brightness from total energy
2. **Outer glow halos** — drawn with `ctx.shadowBlur`, energy-driven blur radius
3. **Rings 7→0** — outside-in so inner rings overlap outer; radius = base + energy * scale; optional wobble via sinusoidal deformation
4. **Center bloom** — bright circle, bass-driven radius and opacity
5. **Particle layer** (optional) — energy-driven sparks emitted from high-energy rings

### 12. Temporal smoothing

Done JS-side: `smooth[t] = alpha * raw[t] + (1 - alpha) * smooth[t-1]`. The smoothing slider controls `alpha`. Different defaults per frequency region could be a future refinement (bass snappier, highs smoother).

### 13. Audio sync strategy

File mode: `<audio>` element plays the WAV. On each `requestAnimationFrame`, read `audio.currentTime`, compute the frame index as `floor(currentTime * fps)`, and look up the corresponding mel data from the precomputed array. Seeking the audio automatically syncs the visual.

### 14. Dev workflow

- Edit `web/` files directly; open `web/index.html` in browser (with local server) for live mic mode
- Run `mel_viz` once to generate an output folder with real data for file-mode testing
- No JS build step — raw `<script type="module">` imports
- `python3 -m http.server` for local development (needed because `file://` blocks ES modules and mic access)
