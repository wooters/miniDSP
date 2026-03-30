## MODIFIED Requirements

### Requirement: Compare script evaluates both VAD systems on LibriVAD test data
The script SHALL load audio files from a LibriVAD Results directory, run both miniDSP VAD and ViT-MFCC inference on each file, and compute F2, precision, recall, and AUC-ROC for both systems against ground-truth labels.

#### Scenario: Full test-clean evaluation
- **WHEN** the user runs `compare_vad.py` with `--librivad-root` pointing to a LibriVAD project containing prepared test-clean data
- **THEN** the script evaluates all 756 files (9 noise types x 6 SNRs x 14 files) and prints overall F2, precision, recall, and AUC-ROC for both miniDSP VAD and ViT-MFCC

#### Scenario: Filtered evaluation
- **WHEN** the user specifies `--noises` and/or `--snrs` to restrict conditions
- **THEN** only matching files are evaluated and metrics (including AUC-ROC) are computed on the subset

### Requirement: Per-condition breakdown output
The script SHALL support a `--breakdown` flag that prints metrics grouped by noise type and by SNR level for both systems.

#### Scenario: Breakdown by noise and SNR
- **WHEN** `--breakdown` is specified
- **THEN** the output includes a per-noise-type table and a per-SNR table showing F2, precision, recall, and AUC-ROC for both systems side by side
