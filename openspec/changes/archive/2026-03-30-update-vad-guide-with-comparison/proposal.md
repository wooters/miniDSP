## Why

The VAD guide (`guides/vad.md`) documents the miniDSP VAD implementation and its Optuna-optimized defaults but stops short of placing the system in context against a published neural baseline. The comparison against ViT-MFCC (small) has been completed (`compare/VAD/README_RESULTS.md`) and reveals important insights — where miniDSP excels, where it falls short, and what the F2/AUC gap means for practitioners choosing between a lightweight C library and a deep-learning model. Including this in the guide helps users make informed deployment decisions without having to find and interpret the raw comparison results separately.

## What Changes

- Add a new section to `guides/vad.md` summarizing the comparison methodology (systems under test, dataset, metrics, frame-rate handling).
- Present the headline results table (F2, precision, recall, AUC-macro, AUC-pooled, wall time) for both miniDSP VAD and ViT-MFCC (small).
- Add per-SNR and per-noise-type breakdowns highlighting where miniDSP is competitive and where the neural model pulls ahead.
- Include the dataset-specific optimization experiment showing that train-clean-100 parameters generalize well (F2 ceiling ~0.935 from threshold/weight tuning alone).
- Add interpretation guidance: what AUC vs F2 reveals about score quality vs. threshold tuning, and the 22× speed advantage of miniDSP.
- Reference the full comparison tooling in `compare/VAD/` for users who want to reproduce or extend the evaluation.

## Capabilities

### New Capabilities

_(none — this change adds content to the existing VAD documentation, not new library capabilities)_

### Modified Capabilities

- `vad-docs`: Adding a "Comparison with ViT-MFCC baseline" section covering methodology, results, per-condition breakdowns, dataset-specific optimization findings, and interpretation guidance.

## Impact

- **Docs**: `guides/vad.md` gains a new section (~100–150 lines). No API or code changes.
- **Spec**: `openspec/specs/vad-docs/spec.md` gains requirements for the new comparison section.
- **No code changes**: The comparison data already exists in `compare/VAD/README_RESULTS.md`; this change curates it into the user-facing guide.
