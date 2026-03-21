## 1. Data Structures and Defaults

- [x] 1.1 Add feature index `#define` constants (`MD_VAD_FEAT_ENERGY` through `MD_VAD_NUM_FEATURES`) to `include/minidsp.h`
- [x] 1.2 Add `MD_vad_params` struct (weights, threshold, onset_frames, hangover_frames, adaptation_rate, band_low_hz, band_high_hz) to `include/minidsp.h`
- [x] 1.3 Add `MD_vad_state` struct (params copy, feat_min/feat_max arrays, onset/hangover counters, current_decision, frames_processed) to `include/minidsp.h`
- [x] 1.4 Add function declarations for `MD_vad_default_params`, `MD_vad_init`, `MD_vad_calibrate`, `MD_vad_process_frame` to `include/minidsp.h`
- [x] 1.5 Create `src/minidsp_vad.c` with `MD_vad_default_params` implementation (equal weights 0.2, threshold 0.5, onset 3, hangover 15, adaptation 0.01, band 300–3400 Hz)
- [x] 1.6 Implement `MD_vad_init` in `src/minidsp_vad.c` (copy params or use defaults, init feat_min to large values, feat_max to small values, zero counters)
- [x] 1.7 Add `test_vad_default_params` test to `tests/test_minidsp.c`

## 2. Feature Extraction Helpers

- [x] 2.1 Implement static `compute_spectral_entropy(psd, num_bins)` — normalize PSD to probability distribution, return -sum(p*log(p)) / log(num_bins)
- [x] 2.2 Implement static `compute_spectral_flatness(psd, num_bins)` — geometric mean / arithmetic mean of PSD bins
- [x] 2.3 Implement static `compute_band_energy_ratio(psd, num_bins, sample_rate, N, band_low_hz, band_high_hz)` — band PSD sum / total PSD sum

## 3. Adaptive Normalization

- [x] 3.1 Implement static `update_normalization(state, raw_features)` — EMA update of feat_min/feat_max with range floor (1e-12)
- [x] 3.2 Implement static `normalize_features(state, raw_features, norm_out)` — map to [0.0, 1.0] using current min/max, clamp result

## 4. Calibration

- [x] 4.1 Implement `MD_vad_calibrate` — compute five raw features, call update_normalization, increment frames_processed, no decision
- [x] 4.2 Add `test_vad_calibrate` test — calibrate on silence frames, verify feat_min/feat_max converge

## 5. Frame Processing and State Machine

- [x] 5.1 Implement `MD_vad_process_frame` — full pipeline: extract features, update normalization, normalize, compute weighted score, apply state machine, write optional outputs, return decision
- [x] 5.2 Add `test_vad_silence` — process silence frames, verify decision stays 0 and score is low
- [x] 5.3 Add `test_vad_speech` — process sine wave in speech band, verify decision becomes 1 after onset
- [x] 5.4 Add `test_vad_hangover` — speech then silence, verify decision holds for hangover_frames then drops
- [x] 5.5 Add `test_vad_onset` — fewer than onset_frames of speech, verify no false trigger
- [x] 5.6 Add `test_vad_features_out` — verify features array is populated, all values in [0.0, 1.0]
- [x] 5.7 Add `test_vad_custom_weights` — set one weight to 1.0, rest to 0.0, verify score tracks that feature

## 6. Build System

- [x] 6.1 Add `src/minidsp_vad.o` to the object list in root `Makefile`

## 7. Doxygen API Documentation

- [x] 7.1 Add full Doxygen doc-comments to all four VAD function declarations in `include/minidsp.h` (LaTeX formulas, `@param`, `@return`, `@code` example, `@see` cross-references)
- [x] 7.2 Add "Voice activity detection (VAD) with adaptive normalization" to the `@brief` feature bullet list in `include/minidsp.h`

## 8. Example Program

- [x] 8.1 Create `examples/vad.c` — synthesize signal with speech/silence segments, init VAD, optional calibration, process frame-by-frame, output CSV, generate interactive HTML (Plotly.js)
- [x] 8.2 Add `//! [vad-init]`, `//! [vad-calibrate]`, `//! [vad-process]`, `//! [vad-custom-weights]` snippet markers
- [x] 8.3 Add `vad` to `EXAMPLES` variable and `plot:` target in `examples/Makefile`
- [x] 8.4 Add `examples/vad`, `examples/vad.csv`, `examples/vad.html` to `.gitignore` and `.dockerignore`

## 9. Interactive HTML Visualization

- [x] 9.1 Build iframe-embeddable HTML visualization showing waveform, feature traces, combined score with threshold, and decision timeline (Plotly.js, 380px iframe height)
- [x] 9.2 Add generated HTML file(s) to `HTML_EXTRA_FILES` in `Doxyfile`

## 10. Tutorial Guide

- [x] 10.1 Create `guides/vad.md` — introduction, feature extraction sections (energy, ZCR, spectral entropy, spectral flatness, band energy ratio) each with LaTeX formula and "Reading the formula in C:" subsection
- [x] 10.2 Add adaptive normalization section with EMA formulas and C reading
- [x] 10.3 Add weighted scoring section with formula and C reading
- [x] 10.4 Add state machine section with transition table
- [x] 10.5 Add visualization iframe embedding with narrative
- [x] 10.6 Add API summary with snippet embedding (`\snippet vad.c vad-init`, etc.)

## 11. Documentation Navigation

- [x] 11.1 Add `\subpage vad` entry to `guides/tutorials.md`
- [x] 11.2 Add VAD to the feature list in `README.md`

## 12. Regenerate llms.txt

- [x] 12.1 Regenerate `llms.txt` and `llms-full.txt` to include new VAD content
