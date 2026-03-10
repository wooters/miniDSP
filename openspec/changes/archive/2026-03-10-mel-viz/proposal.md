## Why

miniDSP has mel-spectrum analysis (`MD_mel_energies`) and STFT, but no visual showcase that demonstrates these in an engaging, non-technical way. An audio visualizer bridges the gap between DSP and something people can *see* and *feel* — and it's a compelling demo of the library's capabilities.

## What Changes

- Add a `tools/` directory as a new home for substantial programs built on miniDSP (distinct from `examples/`, which are API demos)
- Create `tools/mel_viz/` — a mel-spectrum audio visualizer that produces a browser-based radial animation synced to audio playback
- The C program reads a WAV file, computes per-frame mel energies, and assembles an output folder with the mel data + web assets
- The web assets (HTML/JS/CSS) render a Canvas 2D radial visualization with concentric rings, reactive background, bass pulse, and configurable color palettes
- Live microphone mode works without the C program — open `web/index.html` directly and it uses the Web Audio API

## Capabilities

### New Capabilities
- `mel-visualizer`: C program that analyzes a WAV file and produces a folder-based HTML visualization
- `tools-directory`: New `tools/` subdirectory convention for substantial miniDSP-based programs
- `live-mic-mode`: Browser-based mic input using Web Audio API, sharing the same Canvas renderer

### Modified Capabilities
<!-- None — this is entirely additive -->

## Impact

- **Project structure**: New `tools/` directory with its own Makefile conventions
- **Root Makefile**: New `tools` target
- **Dependencies**: No new C dependencies (libsndfile for WAV I/O is already used). Web side is zero-dependency (no npm, no build tools).
- **Build system**: `tools/mel_viz/Makefile` follows the same `config.mk` pattern as `examples/`
