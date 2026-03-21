## Context

miniDSP is a small C library for audio DSP that already provides energy measurement (`MD_energy`), zero-crossing rate (`MD_zero_crossing_rate`), and power spectral density (`MD_power_spectral_density`). These are the standard building blocks for voice activity detection, but the library has no integrated detector. Users who want VAD must wire these primitives together themselves.

The library follows a stateless, caller-owns-memory design: no heap allocations, no global state beyond FFT plan caches, frame-at-a-time processing. The VAD module must follow the same conventions.

## Goals / Non-Goals

**Goals:**
- Provide a frame-at-a-time VAD that composes existing miniDSP primitives into a binary speech/silence decision
- Support adaptive normalization so the detector works across different microphones, gain settings, and noise floors without manual tuning
- Expose all tuning knobs (weights, threshold, onset/hangover frames, adaptation rate, speech band) via a params struct with sensible defaults
- Include onset gating and hangover smoothing for robust decisions
- Maintain the library's zero-heap-allocation, caller-owns-memory contract
- Provide comprehensive documentation: Doxygen API docs, tutorial guide with formulas, example program, interactive visualization

**Non-Goals:**
- Neural/ML-based VAD (e.g., WebRTC VAD, Silero) — this is a classical parametric detector
- Multi-channel or stereo support — single-channel only, consistent with the rest of the library
- Automatic framing — the caller handles windowing and frame extraction, consistent with existing API
- Noise estimation or spectral subtraction — the VAD detects speech presence, not quality
- Real-time audio I/O integration — the example program processes pre-generated or file-based signals

## Decisions

### Five-feature weighted combination over single-feature or ML approaches

Use five complementary features (energy, ZCR, spectral entropy, spectral flatness, band energy ratio) combined via weighted sum. This provides robustness across signal conditions while remaining fully interpretable and tunable.

*Alternative considered*: Single-feature (energy-only) VAD — simpler but fails in moderate noise. ML-based VAD — too heavy for a small C library with no external ML dependencies.

### EMA-based adaptive normalization over fixed thresholds

Track per-feature min/max via exponential moving averages, normalize to [0, 1]. This adapts to different recording conditions automatically.

*Alternative considered*: Fixed normalization ranges — requires manual calibration per deployment. Z-score normalization — needs running variance computation, more state, and unbounded output range.

### Onset gating + hangover over raw thresholding

Require `onset_frames` consecutive above-threshold frames before declaring speech (prevents click false positives). Hold speech decision for `hangover_frames` after score drops (bridges brief dips mid-utterance).

*Alternative considered*: Simple hysteresis (two thresholds) — less intuitive to tune and doesn't address the transient click problem.

### Shared PSD computation for spectral features

Features 2–4 (spectral entropy, flatness, band energy ratio) all derive from the same PSD output. Compute `MD_power_spectral_density` once per frame and pass to all three helpers.

### Static helper functions, no new public API beyond the four entry points

Internal feature extraction helpers (`compute_spectral_entropy`, `compute_spectral_flatness`, `compute_band_energy_ratio`) and normalization helpers (`update_normalization`, `normalize_features`) are `static` in `minidsp_vad.c`. This keeps the public API surface minimal.

### Caller-owned state struct with no heap allocations

`MD_vad_state` is stack-allocatable. All arrays (feat_min, feat_max, etc.) are fixed-size (`MD_VAD_NUM_FEATURES = 5`). No `malloc`/`free` needed. Consistent with the library's existing design.

## Risks / Trade-offs

**[Risk]** PSD computation per frame may be expensive for very short frames on constrained hardware.
→ *Mitigation*: PSD is already cached/planned by `minidsp_spectrum.c`. Users who need lower latency can reduce feature weights to skip spectral features (set weights 2–4 to 0.0) — though the PSD is still computed in the current design. A future optimization could skip PSD when spectral weights are zero.

**[Risk]** EMA adaptation rate is a single parameter controlling all five features — may not suit all features equally.
→ *Mitigation*: Start with a single rate for simplicity. Per-feature rates can be added later as a non-breaking extension to the params struct.

**[Risk]** Fixed feature count (5) baked into struct sizes limits extensibility.
→ *Mitigation*: `MD_VAD_NUM_FEATURES` is a `#define` constant. Adding features requires a version bump but is straightforward.

**[Trade-off]** Band energy ratio uses fixed Hz boundaries (default 300–3400 Hz) — optimized for telephony-band speech. Wideband or non-speech applications may need different bands.
→ *Mitigation*: Band boundaries are caller-tunable params.
