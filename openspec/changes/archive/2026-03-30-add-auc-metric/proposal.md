## Why

The VAD comparison currently reports F2, precision, and recall -- all threshold-dependent metrics. Adding AUC-ROC provides a threshold-independent measure of discriminative ability, showing how well each system separates speech from non-speech across all operating points. This is standard practice in detection system evaluation and will give a more complete picture of each system's capabilities.

## What Changes

- Add an `auc_roc` function that computes the area under the ROC curve from continuous scores and binary targets
- Modify miniDSP VAD evaluation to expose per-frame continuous scores (pre-threshold) alongside binary decisions
- Modify ViT-MFCC evaluation to expose per-frame speech probabilities (pre-argmax) alongside binary decisions
- Add AUC column to overall, per-noise, and per-SNR output tables
- Update README_RESULTS.md with AUC figures

## Capabilities

### New Capabilities
- `auc-metric`: Computation and reporting of AUC-ROC scores for both VAD systems in the comparison framework

### Modified Capabilities
- `vad-comparison`: Add AUC-ROC as an additional metric alongside F2/precision/recall in evaluation output

## Impact

- `compare/VAD/compare_vad.py`: Main changes -- new metric function, modified eval pipelines to return scores, updated print functions
- `compare/VAD/README_RESULTS.md`: Updated with AUC results
- `pyminidsp` (Python bindings): May need to expose raw scores from `VAD.process()` if not already available
- New dependency: `scikit-learn` (for `roc_auc_score`) or hand-rolled trapezoidal AUC
