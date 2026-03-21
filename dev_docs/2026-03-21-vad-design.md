# Voice Activity Detection (VAD) — Design Overview

## Motivation

miniDSP already provides the building blocks for voice activity detection — energy, zero-crossing rate, spectral analysis — but no integrated detector. Adding a parametric VAD gives users a practical, real-time-capable utility and a clear example of how these primitives compose into a decision system.

## Design Summary

A frame-at-a-time VAD that combines five normalized audio features into a weighted score, applies a threshold, and smooths the binary decision with onset gating and hangover.

### Features

| # | Feature           | Measures                        | Existing API used            |
|---|-------------------|---------------------------------|------------------------------|
| 0 | Energy            | Frame loudness                  | `MD_energy`                  |
| 1 | Zero-crossing rate| Temporal structure               | `MD_zero_crossing_rate`      |
| 2 | Spectral entropy  | Spectral randomness             | Derived from `MD_power_spectral_density` |
| 3 | Spectral flatness | Tonal vs. noisy spectrum        | Derived from `MD_power_spectral_density` |
| 4 | Band energy ratio | Speech-band concentration       | Derived from `MD_power_spectral_density` |

Features 2–4 share one PSD computation per frame.

### Adaptive Normalization

Raw feature values are normalized to [0.0, 1.0] using exponential-moving-average estimates of each feature's min and max. This lets the detector adapt to different microphones, gain settings, and noise floors without manual tuning. The adaptation rate is a caller-tunable parameter.

### State Machine

The detector uses an onset + hangover mechanism:

- **Onset**: the combined score must exceed the threshold for `onset_frames` consecutive frames before the decision flips to speech. This prevents transient clicks from triggering false positives.
- **Hangover**: after the score drops below threshold, the speech decision holds for `hangover_frames` additional frames. This bridges brief dips mid-utterance.

### API Shape

- `MD_vad_default_params` — populate a params struct with sensible defaults.
- `MD_vad_init` — initialize state from params (or NULL for defaults).
- `MD_vad_calibrate` — optional: feed known-silence frames to seed the adaptive normalization.
- `MD_vad_process_frame` — process one frame, return binary decision, optionally output the combined score and per-feature breakdown.

All state lives in a caller-owned `MD_vad_state` struct (no heap allocations). The caller handles framing — consistent with the rest of the library.

### File Layout

- Structs and declarations: `include/minidsp.h`
- Implementation: `src/minidsp_vad.c`
- Tests: `tests/test_minidsp.c`
