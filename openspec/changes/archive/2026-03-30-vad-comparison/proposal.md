## Why

After optimizing the miniDSP VAD hyperparameters on LibriVAD train-clean-100 (F2: 0.837 → 0.933), we need to understand how the system performs against a state-of-the-art neural baseline. The LibriVAD project provides a Vision Transformer (ViT-MFCC) baseline — a ~3.5M parameter deep learning model — making it an ideal reference point. Comparing a lightweight, real-time C library against a large trained model on the same test data quantifies the gap (or lack thereof) and guides future work.

## What Changes

- Add a `compare/VAD/` directory with a standalone Python script that evaluates both miniDSP VAD and ViT-MFCC on the LibriVAD test-clean split
- The script downloads the pre-trained ViT-MFCC (small) model from HuggingFace automatically
- Both systems are evaluated on the same 756 test files (9 noise types × 6 SNRs × 14 files) using F2, precision, and recall
- Results are reported as overall metrics and per-condition breakdowns (by noise type and SNR)
- A README documents setup, usage, and how to interpret results

## Capabilities

### New Capabilities
- `vad-comparison`: Standalone comparison framework that runs miniDSP VAD and ViT-MFCC side-by-side on LibriVAD test data, computing common metrics (F2, precision, recall) for both systems

### Modified Capabilities
<!-- None -->

## Impact

- **New files**: `compare/VAD/compare_vad.py`, `compare/VAD/pyproject.toml`, `compare/VAD/README.md`
- **Dependencies**: New Python project with torch, huggingface_hub, einops, scipy, numpy, scikit-learn, soundfile, pyminidsp
- **External data**: Requires LibriVAD repo (sibling directory) with test-clean data prepared, and imports ViT code from `LibriVAD/Eval-ViT-MFCC/`
- **No changes** to the miniDSP C library, headers, build system, or existing optimize/ scripts
