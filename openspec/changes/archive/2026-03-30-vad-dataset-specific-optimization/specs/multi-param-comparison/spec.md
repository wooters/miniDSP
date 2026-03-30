## ADDED Requirements

### Requirement: Load VAD parameters from JSON files

The comparison script SHALL support loading optimized VAD parameters from JSON files produced by `optimize_vad.py`. The JSON format contains `weights_normalized` (array of 5 floats) and `params` (object with `threshold`, `onset_frames`, `hangover_frames`, `adaptation_rate`, `band_low_hz`, `band_high_hz`).

#### Scenario: Valid parameter file loaded
- **WHEN** a valid `best_params.json` file path is provided
- **THEN** the script SHALL construct a `VAD()` instance using the weights and parameters from that file

#### Scenario: Missing parameter file
- **WHEN** a parameter file path does not exist
- **THEN** the script SHALL exit with a clear error message naming the missing file

### Requirement: Multiple named parameter sets via CLI

The comparison script SHALL accept a repeatable `--params` argument in `label=path` format, where `label` is the display name and `path` is a JSON parameter file.

#### Scenario: Multiple parameter sets specified
- **WHEN** the user provides `--params "train-clean-100=path1.json" --params "dev-clean=path2.json"`
- **THEN** the script SHALL evaluate each parameter set independently on the test data and include all results in the output

#### Scenario: No params specified (backward compatibility)
- **WHEN** no `--params` arguments are provided
- **THEN** the script SHALL evaluate using library default parameters labeled "miniDSP VAD", preserving current behavior

### Requirement: Side-by-side output tables

The overall and breakdown output tables SHALL display metrics for all miniDSP parameter sets alongside the ViT-MFCC baseline.

#### Scenario: Overall table with three parameter sets
- **WHEN** three miniDSP parameter sets are evaluated
- **THEN** the overall table SHALL show one row per miniDSP configuration plus one row for ViT-MFCC, each with F-beta, precision, recall, AUC(macro), AUC(pooled), and timing

#### Scenario: Breakdown tables with multiple parameter sets
- **WHEN** `--breakdown` is enabled with multiple parameter sets
- **THEN** per-noise and per-SNR tables SHALL include one column per miniDSP configuration and one column for ViT-MFCC

### Requirement: Independent evaluation per parameter set

Each miniDSP parameter set SHALL be evaluated independently with a fresh `VAD()` instance per file, ensuring no state leaks between configurations.

#### Scenario: Parameter sets produce different results
- **WHEN** two parameter sets with different thresholds are evaluated on the same data
- **THEN** their F-beta scores SHALL differ (unless parameters happen to produce identical decisions)

#### Scenario: Each file gets a fresh VAD state
- **WHEN** a parameter set is evaluated across multiple files
- **THEN** each file SHALL be processed with a newly initialized `VAD()` instance (no state carried between files)
