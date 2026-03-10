# mel_viz

A mel-spectrum audio visualizer that generates browser-based radial animations synced to audio playback.

## Build

From the repository root:

```bash
make tools
```

Or directly:

```bash
make -C tools/mel_viz
```

## Usage

### File mode (analyze a WAV file)

```bash
./tools/mel_viz/mel_viz song.wav -o output/
open output/index.html
```

The output folder contains a self-contained web visualization — just open `index.html` in any browser.

### Live mic mode

Open `tools/mel_viz/web/index.html` directly via a local server:

```bash
cd tools/mel_viz
python3 -m http.server 8000
# open http://localhost:8000/web/
```

Grant microphone access when prompted. No WAV file or C program needed.

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-o <dir>` | `mel_viz_out/` | Output directory |
| `--mels <n>` | 24 | Number of mel bands |
| `--fft-size <n>` | 2048 | FFT window size |
| `--fps <n>` | 30 | Frames per second |
| `--min-freq <hz>` | 40 | Low frequency bound |
| `--max-freq <hz>` | 16000 | High frequency bound |
| `--width <px>` | 1080 | Canvas width |
| `--height <px>` | 1080 | Canvas height |

## Visual Controls

The browser UI provides real-time knobs (no recomputation needed):

- **Palette**: Plasma, Ocean, Fire, Neon
- **Rings**: 8, 12, or 24 visual ring groups
- **Smoothing**: Temporal smoothing (snappy to smooth)
- **Bass**: Bass pulse sensitivity (global zoom on kicks)
- **Wobble**: Ring deformation amount
- **Glow**: Neon glow intensity

## Architecture

```
WAV file --> mel_viz (C) --> data.js + web assets --> Browser (Canvas 2D)
                                                         |
Microphone --> Web Audio API --> JS mel weighting --------+
                                                         |
                                                    Canvas renderer
                                                    (radial rings,
                                                     reactive background,
                                                     bass pulse, glow)
```

The C program handles DSP (reading audio, computing mel energies via `MD_mel_energies()`). The browser handles rendering via Canvas 2D. File mode and mic mode share the same renderer — only the data source differs.
