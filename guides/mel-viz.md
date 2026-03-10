# mel_viz -- Mel-Spectrum Audio Visualizer {#mel-viz}

**mel_viz** is a tool that generates browser-based radial animations driven by mel-spectrum analysis.
Feed it a WAV file and it produces an interactive HTML visualization synced to audio playback.
It also supports live microphone input directly in the browser.

## Build

From the repository root:

```sh
make tools
```

## File mode

Analyze a WAV file and generate a visualization folder:

```sh
./tools/mel_viz/mel_viz samples/punchy_slap_bass_30s.wav -o /tmp/viz
cd /tmp/viz && python3 -m http.server 8000
# open http://localhost:8000
```

The `samples/` directory at the repository root contains audio files
ready to use with mel_viz (e.g., `punchy_slap_bass_30s.wav`).

The output folder is self-contained -- open `index.html` in any browser to
see the visualization with audio playback, palette selection, and real-time
visual knobs.

## Live mic mode

The web renderer works standalone without the C program.
Serve the `web/` directory and click **Mic**:

```sh
cd tools/mel_viz && python3 -m http.server 8000
# open http://localhost:8000/web/
```

Grant microphone access when prompted.
The visualization responds to live audio using the Web Audio API.

## Architecture

```
File mode:
  WAV file --> mel_viz (C) --> data.js + web assets --> Browser (Canvas 2D)
                 |                                          |
                 | MD_mel_energies()                   audio player sync
                 | (per-frame mel analysis)            + radial renderer

Mic mode:
  Microphone --> Web Audio API --> JS mel weighting --> Same renderer
                 getUserMedia()    AnalyserNode
```

The C program reads audio via libsndfile, computes per-frame mel energies
using `MD_mel_energies()`, and writes the data as a JavaScript file.
The browser renders concentric rings whose size, color, glow, and wobble
are driven by the mel band energies.
A bass pulse effect scales the entire canvas on kick/bass hits.

**Visual knobs** (adjustable in real time, no recomputation):
palette, smoothing, bass sensitivity, ring wobble, glow intensity, and ring count.

For the full list of CLI flags and visual controls, see the
[mel_viz README](https://github.com/wooters/miniDSP/blob/main/tools/mel_viz/README.md).
