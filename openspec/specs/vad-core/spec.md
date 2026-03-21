## ADDED Requirements

### Requirement: VAD feature constants
The library SHALL define `MD_VAD_FEAT_ENERGY` (0), `MD_VAD_FEAT_ZCR` (1), `MD_VAD_FEAT_SPECTRAL_ENTROPY` (2), `MD_VAD_FEAT_SPECTRAL_FLATNESS` (3), `MD_VAD_FEAT_BAND_ENERGY_RATIO` (4), and `MD_VAD_NUM_FEATURES` (5) as preprocessor constants in `include/minidsp.h`.

#### Scenario: Feature index constants are defined
- **WHEN** a program includes `minidsp.h`
- **THEN** all six `MD_VAD_FEAT_*` and `MD_VAD_NUM_FEATURES` constants SHALL be available as integer preprocessor defines

### Requirement: VAD params struct with sensible defaults
The library SHALL provide an `MD_vad_params` struct containing: `weights[MD_VAD_NUM_FEATURES]`, `threshold`, `onset_frames`, `hangover_frames`, `adaptation_rate`, `band_low_hz`, `band_high_hz`. The function `MD_vad_default_params` SHALL populate a params struct with: equal weights (0.2 each), threshold 0.5, onset 3 frames, hangover 15 frames, adaptation rate 0.01, band 300–3400 Hz.

#### Scenario: Default params are populated correctly
- **WHEN** `MD_vad_default_params` is called with a valid params pointer
- **THEN** all fields SHALL be set to the documented default values

#### Scenario: Null params pointer
- **WHEN** `MD_vad_default_params` is called with a NULL pointer
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER` and the function SHALL return without crashing

### Requirement: VAD state initialization
`MD_vad_init` SHALL initialize a caller-owned `MD_vad_state` struct from a params struct. If the params pointer is NULL, it SHALL use default params. After initialization, feat_min SHALL be set to large values, feat_max to small values, counters and decision to zero, and frames_processed to zero.

#### Scenario: Init with explicit params
- **WHEN** `MD_vad_init` is called with a valid state pointer and a non-NULL params pointer
- **THEN** the state SHALL copy the provided params and initialize all internal fields

#### Scenario: Init with NULL params (use defaults)
- **WHEN** `MD_vad_init` is called with a valid state pointer and NULL params
- **THEN** the state SHALL use default params from `MD_vad_default_params`

#### Scenario: Null state pointer
- **WHEN** `MD_vad_init` is called with a NULL state pointer
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER`

### Requirement: Adaptive normalization via EMA
The VAD SHALL maintain per-feature min and max estimates updated via exponential moving average on each frame. Features SHALL be normalized to [0.0, 1.0] using the current min/max range. A minimum range floor (e.g., 1e-12) SHALL be enforced to prevent division by zero.

#### Scenario: Normalization adapts over time
- **WHEN** frames with varying feature values are processed
- **THEN** feat_min and feat_max SHALL converge toward the observed range and normalized features SHALL fall within [0.0, 1.0]

#### Scenario: Zero-range protection
- **WHEN** all frames have identical feature values (min == max)
- **THEN** normalization SHALL NOT produce NaN or Inf — the range floor SHALL prevent division by zero

### Requirement: VAD calibration on known silence
`MD_vad_calibrate` SHALL compute all five raw features on an input frame, update adaptive normalization min/max estimates, increment frames_processed, and return without producing a decision or modifying the state machine.

#### Scenario: Calibrate seeds normalization
- **WHEN** several known-silence frames are passed to `MD_vad_calibrate`
- **THEN** feat_min and feat_max SHALL converge to values representative of silence

#### Scenario: Calibrate does not affect decision state
- **WHEN** frames are passed only to `MD_vad_calibrate` (never to `MD_vad_process_frame`)
- **THEN** `current_decision` SHALL remain 0 and onset/hangover counters SHALL remain at their initial values

#### Scenario: Null pointer arguments
- **WHEN** `MD_vad_calibrate` is called with NULL state or NULL signal pointer
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER`

### Requirement: VAD frame processing with weighted scoring
`MD_vad_process_frame` SHALL compute five raw features (energy, ZCR, spectral entropy, spectral flatness, band energy ratio), update adaptive normalization, normalize features to [0.0, 1.0], compute a weighted score as `sum(weights[i] * norm[i])`, apply the onset/hangover state machine, and return a binary decision (0 or 1).

#### Scenario: Silence produces decision 0
- **WHEN** silence frames (all zeros) are processed
- **THEN** the decision SHALL remain 0 and the combined score SHALL be low

#### Scenario: Speech-like signal produces decision 1
- **WHEN** sine wave frames in the speech band (300–3400 Hz) are processed for enough frames to exceed onset
- **THEN** the decision SHALL become 1 after the onset period

#### Scenario: Optional score output
- **WHEN** `MD_vad_process_frame` is called with a non-NULL `score_out` pointer
- **THEN** the combined weighted score SHALL be written to `*score_out`

#### Scenario: Optional features output
- **WHEN** `MD_vad_process_frame` is called with a non-NULL `features_out` array
- **THEN** all five normalized feature values SHALL be written, each in [0.0, 1.0]

#### Scenario: NULL score_out and features_out are allowed
- **WHEN** `MD_vad_process_frame` is called with NULL `score_out` and/or NULL `features_out`
- **THEN** the function SHALL skip writing those outputs without error

#### Scenario: Null pointer for required arguments
- **WHEN** `MD_vad_process_frame` is called with NULL state or NULL signal pointer
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER` and a safe default (0) SHALL be returned

### Requirement: Onset gating prevents false triggers
The state machine SHALL require `onset_frames` consecutive above-threshold frames before transitioning from silence to speech. If the score drops below threshold before reaching `onset_frames`, the onset counter SHALL reset.

#### Scenario: Transient click does not trigger speech
- **WHEN** fewer than `onset_frames` consecutive above-threshold frames are followed by below-threshold frames
- **THEN** the decision SHALL remain 0

#### Scenario: Sustained signal triggers after onset period
- **WHEN** `onset_frames` or more consecutive above-threshold frames are processed
- **THEN** the decision SHALL transition to 1

### Requirement: Hangover prevents premature silence
After the score drops below threshold during speech, the speech decision SHALL hold for `hangover_frames` additional frames before transitioning to silence.

#### Scenario: Brief dip mid-utterance is bridged
- **WHEN** speech frames are followed by a few below-threshold frames (fewer than `hangover_frames`) then above-threshold frames resume
- **THEN** the decision SHALL remain 1 throughout

#### Scenario: Sustained silence after hangover causes transition
- **WHEN** speech frames are followed by `hangover_frames` or more below-threshold frames
- **THEN** the decision SHALL transition to 0

### Requirement: Custom weights control feature influence
The caller SHALL be able to set arbitrary weights in `MD_vad_params.weights`. The weighted score SHALL be computed as the dot product of weights and normalized features.

#### Scenario: Single-feature mode
- **WHEN** one weight is set to 1.0 and all others to 0.0
- **THEN** the combined score SHALL track only the selected feature's normalized value

### Requirement: Shared PSD computation
Features 2–4 (spectral entropy, spectral flatness, band energy ratio) SHALL share a single `MD_power_spectral_density` computation per frame.

#### Scenario: PSD computed once per frame
- **WHEN** `MD_vad_process_frame` is called
- **THEN** `MD_power_spectral_density` SHALL be called exactly once, and the result SHALL be passed to all three spectral feature helpers

### Requirement: No heap allocations
All VAD state SHALL reside in the caller-owned `MD_vad_state` struct. The VAD module SHALL NOT call `malloc`, `calloc`, `realloc`, or `free`.

#### Scenario: Stack-allocated state works
- **WHEN** `MD_vad_state` is declared on the stack and initialized via `MD_vad_init`
- **THEN** all VAD operations SHALL function correctly without any heap allocation

### Requirement: Build system integration
`src/minidsp_vad.o` SHALL be added to the object list for `libminidsp.a` in the root Makefile.

#### Scenario: Library builds with VAD module
- **WHEN** `make` is run at the repo root
- **THEN** `src/minidsp_vad.o` SHALL be compiled and linked into `libminidsp.a`
