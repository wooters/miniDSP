## 1. Dependencies

- [x] 1.1 Add `scikit-learn` to `compare/VAD/pyproject.toml`

## 2. Score Plumbing

- [x] 2.1 Change `eval_minidsp` to capture and return scores from `VAD.process()` — return `list[tuple[preds, targets, scores]]`
- [x] 2.2 Change `vit_infer_file` to compute `softmax(logits, dim=1)[:, 1]` and return speech probabilities alongside predictions
- [x] 2.3 Change `eval_vit` to thread scores through — return `list[tuple[preds, targets, scores]]`

## 3. AUC Computation

- [x] 3.1 Add `compute_auc` function that wraps `sklearn.metrics.roc_auc_score` with single-class edge-case handling (return `None` when only one class present)
- [x] 3.2 Update `aggregate_metrics` to accept optional scores array and compute AUC
- [x] 3.3 Update `per_condition_metrics` to thread scores through to `aggregate_metrics`

## 4. Output Tables

- [x] 4.1 Update `print_overall` to display AUC column
- [x] 4.2 Update `print_breakdown` to display AUC columns in per-noise and per-SNR tables

## 5. Main Function

- [x] 5.1 Update `main()` to pass scores through the pipeline to print functions

## 6. Run Comparison and Update Docs

- [x] 6.1 Re-run the full VAD comparison (`uv run python compare_vad.py --librivad-root ... --breakdown`) and capture output
- [x] 6.2 Update `README_RESULTS.md`: add AUC column to Overall Results table, Per SNR table, Per Noise Type table, and update Key Takeaways with AUC observations
- [x] 6.3 Update `README.md`: add AUC to the "Metric definitions" section and mention AUC in the "Interpreting results" / "ViT threshold" sections
