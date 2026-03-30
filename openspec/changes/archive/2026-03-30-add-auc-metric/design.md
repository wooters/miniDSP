## Context

The VAD comparison (`compare/VAD/compare_vad.py`) currently evaluates miniDSP VAD and ViT-MFCC using only threshold-dependent metrics (F2, precision, recall). Both systems already produce continuous scores internally:

- **miniDSP VAD**: `VAD.process()` returns `(decisions, scores, features)` where `scores` is a float64 array of weighted combined scores in [0.0, 1.0]. The comparison currently discards `scores` (line 389: `decisions, _, _ = vad.process(...)`).
- **ViT-MFCC**: The KWT model produces logits, from which `torch.softmax(..., dim=1)[:, 1]` gives speech probabilities in [0.0, 1.0]. The comparison currently applies `torch.argmax` and discards probabilities.

## Goals / Non-Goals

**Goals:**
- Add AUC-ROC as a threshold-independent metric to the overall, per-noise, and per-SNR output tables
- Expose continuous scores from both eval pipelines without changing existing binary decision outputs
- Keep the change self-contained to `compare_vad.py` and results documentation

**Non-Goals:**
- Adding other threshold-independent metrics (e.g., AP, DET curves, EER)
- Changing the miniDSP C library or pyminidsp bindings (scores are already exposed)
- Adding ROC curve plotting

## Decisions

### 1. Use sklearn's `roc_auc_score` vs hand-rolled trapezoidal AUC

**Decision**: Use `sklearn.metrics.roc_auc_score`.

**Rationale**: scikit-learn is the standard and handles edge cases (single-class inputs, ties). The compare project already has heavyweight dependencies (torch, einops, huggingface_hub), so adding sklearn is negligible. A hand-rolled version saves one dependency but risks subtle correctness issues.

**Alternative considered**: `np.trapz` on manually sorted FPR/TPR pairs. Rejected because edge-case handling adds complexity without meaningful benefit.

### 2. Return type change for eval functions

**Decision**: Change `eval_minidsp` and `eval_vit` return types from `list[tuple[predictions, targets]]` to `list[tuple[predictions, targets, scores]]` — adding continuous scores as a third element.

**Rationale**: This is the minimal change that threads scores through the existing aggregation pipeline. Callers that only need binary decisions can ignore the third element. Downstream functions (`aggregate_metrics`, `per_condition_metrics`) gain an optional `compute_auc` parameter.

**Alternative considered**: Separate `eval_minidsp_scores` functions. Rejected because it would duplicate the evaluation loop.

### 3. AUC for miniDSP scores

**Decision**: Use the raw weighted score from `VAD.process()` directly as the AUC score input.

**Rationale**: The score is already a continuous value in [0, 1] that represents the pre-threshold detection confidence. It is the natural operating-point variable for ROC analysis.

### 4. AUC for ViT-MFCC scores

**Decision**: Apply `torch.softmax(logits, dim=1)[:, 1]` to get per-frame speech probability and use that as the AUC score input.

**Rationale**: Softmax probability of the speech class is the standard way to get continuous scores from a classification model. This is what `argmax` thresholds at 0.5.

## Risks / Trade-offs

- **[Single-class edge case]** If a subset (e.g., one noise@SNR bucket) has all-speech or all-silence targets, `roc_auc_score` raises `ValueError`. **Mitigation**: Catch and report `N/A` for that subset.
- **[Minimal performance impact]** Softmax computation adds negligible overhead to ViT inference. Storing score arrays alongside decisions roughly doubles memory for results, but total frame counts (~750k) make this trivial.
- **[New dependency]** Adding `scikit-learn` to `pyproject.toml`. **Mitigation**: It's a standard, well-maintained package already in the PyPI ecosystem. The compare project is not a lightweight tool — it already requires torch.
