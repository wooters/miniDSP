## Why

miniDSP provides the building blocks for voice activity detection — energy, zero-crossing rate, spectral analysis — but no integrated detector. Adding a parametric VAD gives users a practical, real-time-capable utility and a clear example of how these primitives compose into a decision system. VAD is a fundamental preprocessing step for streaming audio, ASR pipelines, and noise reduction.

## What Changes

- New `MD_vad_*` public API: `MD_vad_default_params`, `MD_vad_init`, `MD_vad_calibrate`, `MD_vad_process_frame`
- New `MD_vad_params` and `MD_vad_state` structs for caller-owned, heap-free state management
- Five normalized audio features (energy, ZCR, spectral entropy, spectral flatness, band energy ratio) combined via weighted scoring
- Adaptive normalization using EMA-tracked min/max per feature
- Onset gating and hangover smoothing state machine for robust binary decisions
- New implementation file `src/minidsp_vad.c` and build system integration
- Example program `examples/vad.c` with CSV output and interactive HTML visualization
- Tutorial guide `guides/vad.md` with formulas, code walkthroughs, and embedded visualizations
- Interactive HTML visualization asset for documentation

## Capabilities

### New Capabilities
- `vad-core`: Frame-at-a-time VAD engine — feature extraction, adaptive normalization, weighted scoring, onset/hangover state machine, and the four public API functions
- `vad-docs`: Tutorial guide, example program, HTML visualization, Doxygen integration, and README/navigation updates

### Modified Capabilities

(none)

## Impact

- **New files**: `src/minidsp_vad.c`, `examples/vad.c`, `guides/vad.md`
- **Modified files**: `include/minidsp.h` (structs, declarations, Doxygen), `Makefile` (new object), `examples/Makefile`, `tests/test_minidsp.c`, `guides/tutorials.md`, `README.md`, `Doxyfile`, `.gitignore`, `.dockerignore`
- **Dependencies**: No new external dependencies — uses existing FFTW3 (via `MD_power_spectral_density`), `MD_energy`, `MD_zero_crossing_rate`
- **API**: Additive only — no breaking changes
