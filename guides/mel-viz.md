# mel_viz -- Mel-Spectrum Audio Visualizer {#mel-viz}

**mel_viz** turns any WAV file into a browser-based radial animation driven by
mel-spectrum analysis.  Concentric rings pulse, wobble, and glow in response to
the audio — bass hits zoom the whole canvas.

\htmlonly
<video controls loop playsinline width="480"
       style="display:block; margin:1em auto; border-radius:8px;">
  <source src="mel-viz-demo.mp4" type="video/mp4">
  Your browser does not support the video tag.
</video>
\endhtmlonly

\image html mel-viz-screenshot.png "mel_viz running on a bass-heavy track (Plasma palette)" width=600px

## Visualize your own audio

Follow these steps to create a visualization from any WAV file.

**1. Build mel_viz**

From the repository root:

```sh
make tools
```

This compiles `tools/mel_viz/mel_viz` (requires FFTW3 and libsndfile).

**2. Run mel_viz on your WAV file**

```sh
./tools/mel_viz/mel_viz your_song.wav -o /tmp/viz
```

Replace `your_song.wav` with the path to your audio file.  The tool reads the
WAV, computes per-frame mel energies, and writes a self-contained visualization
folder to the output directory.

Common options:

| Flag | Default | What it does |
|------|---------|--------------|
| `-o <dir>` | `mel_viz_out/` | Output directory |
| `--mels <n>` | 24 | Number of mel bands |
| `--fft-size <n>` | 2048 | FFT window size |
| `--fps <n>` | 30 | Frames per second |
| `--min-freq <hz>` | 40 | Low frequency bound |
| `--max-freq <hz>` | 16000 | High frequency bound |

**3. Open in a browser**

The output folder needs to be served over HTTP (browsers block local `file://`
audio playback).  The easiest way:

```sh
cd /tmp/viz
python3 -m http.server 8000
```

Then open <http://localhost:8000> in your browser and press play.

**4. Tweak the visual controls**

The side panel on the left has real-time knobs — no recomputation needed:

- **Palette** — Plasma, Ocean, Fire, or Neon color themes
- **Rings** — 8, 12, or 24 concentric ring groups
- **Smoothing** — Temporal smoothing from snappy to flowing
- **Bass** — How much kick/bass hits zoom the entire canvas
- **Wobble** — Ring edge deformation amount
- **Glow** — Neon bloom intensity around each ring

Experiment with different palettes and cranking up the glow and wobble for a
more psychedelic look, or dial them down for a clean, minimal animation.

**5. Export as video**

Want to share your visualization?  Click **Export MP4** in the side panel.
The browser renders every frame offline and muxes it into an MP4 file
(Chrome or Edge required — uses the WebCodecs API).

## Live mic mode

The web renderer also works standalone without the C program — it captures
microphone audio and visualizes it in real time.

```sh
cd tools/mel_viz
python3 -m http.server 8000
# open http://localhost:8000/web/
```

Click **Mic**, grant microphone access, and the visualization responds to
whatever you play or say.

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

For the full list of CLI flags, see the
[mel_viz README](https://github.com/wooters/miniDSP/blob/main/tools/mel_viz/README.md).
