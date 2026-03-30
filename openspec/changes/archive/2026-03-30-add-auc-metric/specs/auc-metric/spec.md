## ADDED Requirements

### Requirement: AUC-ROC computed from continuous scores and binary targets
The comparison script SHALL compute the area under the receiver operating characteristic curve (AUC-ROC) using `sklearn.metrics.roc_auc_score` from per-frame continuous scores and binary ground-truth targets.

#### Scenario: AUC from miniDSP weighted scores
- **WHEN** miniDSP VAD is evaluated on a set of files
- **THEN** the per-frame weighted scores (float64 in [0.0, 1.0]) from `VAD.process()` are collected alongside binary decisions, and `roc_auc_score(targets, scores)` is computed on the concatenated arrays

#### Scenario: AUC from ViT-MFCC speech probabilities
- **WHEN** ViT-MFCC is evaluated on a set of files
- **THEN** per-frame speech probabilities are obtained via `softmax(logits, dim=1)[:, 1]` and collected alongside binary decisions, and `roc_auc_score(targets, scores)` is computed on the concatenated arrays

#### Scenario: Single-class subset
- **WHEN** a subset of evaluation data contains only one class (all speech or all silence)
- **THEN** AUC SHALL be reported as `N/A` instead of raising an error

### Requirement: miniDSP evaluation exposes continuous scores
The `eval_minidsp` function SHALL return per-frame continuous scores in addition to binary decisions and targets.

#### Scenario: Scores captured from VAD.process
- **WHEN** `eval_minidsp` processes a list of audio files
- **THEN** it returns a list of `(predictions, targets, scores)` tuples where `scores` is the float64 array from `VAD.process()`

### Requirement: ViT-MFCC evaluation exposes continuous scores
The `vit_infer_file` and `eval_vit` functions SHALL return per-frame speech probabilities in addition to binary decisions and targets.

#### Scenario: Probabilities captured from softmax
- **WHEN** `vit_infer_file` processes an audio file
- **THEN** it applies `softmax(logits, dim=1)[:, 1]` to get per-frame speech probability and returns `(predictions, targets, scores)`

### Requirement: AUC displayed in output tables
The comparison output tables SHALL include an AUC column alongside F-beta, precision, and recall.

#### Scenario: Overall table includes AUC
- **WHEN** the overall comparison table is printed
- **THEN** it includes an `AUC` column showing the AUC-ROC for each system

#### Scenario: Breakdown tables include AUC
- **WHEN** `--breakdown` is specified
- **THEN** per-noise and per-SNR tables include miniDSP AUC and ViT AUC columns
