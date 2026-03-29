## MODIFIED Requirements

### Requirement: VAD params struct with sensible defaults
The library SHALL provide an `MD_vad_params` struct containing: `weights[MD_VAD_NUM_FEATURES]`, `threshold`, `onset_frames`, `hangover_frames`, `adaptation_rate`, `band_low_hz`, `band_high_hz`. The function `MD_vad_default_params` SHALL populate a params struct with the following Optuna-optimized values (300 trials, F2 metric, LibriVAD train-clean-100):

| Parameter | Value |
|---|---|
| `weights[0]` (energy) | 0.723068 |
| `weights[1]` (zcr) | 0.063948 |
| `weights[2]` (spectral_entropy) | 0.005964 |
| `weights[3]` (spectral_flatness) | 0.048865 |
| `weights[4]` (band_energy_ratio) | 0.158156 |
| `threshold` | 0.245332 |
| `onset_frames` | 1 |
| `hangover_frames` | 22 |
| `adaptation_rate` | 0.012755 |
| `band_low_hz` | 126.4 |
| `band_high_hz` | 2899.3 |

#### Scenario: Default params are populated correctly
- **WHEN** `MD_vad_default_params` is called with a valid params pointer
- **THEN** all fields SHALL be set to the optimized default values listed above

#### Scenario: Null params pointer
- **WHEN** `MD_vad_default_params` is called with a NULL pointer
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER` and the function SHALL return without crashing
