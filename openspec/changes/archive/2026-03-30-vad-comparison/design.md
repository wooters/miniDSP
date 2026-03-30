## Context

The miniDSP VAD is a lightweight, real-time voice activity detector implemented in C (~300 lines), using 5 hand-crafted features with a weighted threshold and onset/hangover state machine. Its hyperparameters were recently optimized via Optuna on LibriVAD train-clean-100, achieving F2=0.933.

The ViT-MFCC baseline from the LibriVAD project is a Vision Transformer (KWT architecture) with ~3.5M parameters, operating on 39-dimensional MFCC features. Pre-trained model checkpoints are available on HuggingFace.

Both systems can be evaluated on LibriVAD test-clean (756 files, 9 noise types × 6 SNR levels), but they currently report different metrics. A fair comparison requires running both on the same data with the same metrics.

## Goals / Non-Goals

**Goals:**
- Produce F2, precision, and recall numbers for both systems on LibriVAD test-clean
- Support per-condition breakdowns (by noise type and by SNR)
- Auto-download the ViT-MFCC checkpoint from HuggingFace
- Self-contained setup via `uv` with a clear README

**Non-Goals:**
- Training or fine-tuning the ViT model
- Optimizing the ViT threshold (use the designed 0.5 operating point)
- Adding AUC/EER metrics (can be added later)
- Modifying the miniDSP C library or its defaults
- Supporting datasets other than LibriSpeechConcat

## Decisions

### 1. Import ViT code from sibling LibriVAD repo

The KWT model class, MFCC extraction, and sequence utilities live in `../LibriVAD/Eval-ViT-MFCC/`. The script will add this path to `sys.path` at runtime using the `--librivad-root` argument.

**Why not vendor the code:** The ViT implementation is the reference — using it directly ensures credibility and avoids maintaining a copy. The LibriVAD repo is already required for the test data anyway.

**Why not install as a package:** The ViT code has no `setup.py` or packaging; it's a collection of scripts.

### 2. Each system evaluated at its native frame rate

miniDSP operates at 20ms frames (50 fps); ViT-MFCC at 10ms frames (100 fps). Ground-truth labels are downsampled from sample-level to each system's frame rate independently via majority vote.

**Why not force a common frame rate:** Forcing both to 20ms would artificially handicap the ViT's temporal resolution, which is an inherent design advantage. Evaluating at native rates reflects real-world performance.

### 3. ViT threshold fixed at 0.5

The ViT outputs 2-class softmax probabilities. The decision boundary is `argmax`, equivalent to threshold=0.5.

**Why not optimize the ViT threshold:** The model was trained with cross-entropy loss where 0.5 is the designed operating point. The miniDSP threshold *was* optimized, but that's an inherent advantage of the miniDSP system's simplicity — its threshold is a tunable knob. Using 0.5 for ViT is the fair, standard choice.

### 4. Apple Silicon MPS acceleration

PyTorch on macOS supports Metal Performance Shaders for GPU acceleration. The script will prefer MPS > CUDA > CPU automatically.

### 5. Directory structure mirrors optimize/VAD/

Place the comparison in `compare/VAD/` with its own `pyproject.toml` (uv project) and README, following the same pattern as `optimize/VAD/`.

### 6. Metric computation reused from optimize_vad.py

The F-beta computation (TP/FP/FN accumulation → precision/recall/F2) is copied from `optimize_vad.py` rather than factored into a shared library, keeping both scripts self-contained.

## Risks / Trade-offs

**ViT import fragility** → The script documents the required LibriVAD repo location and fails with a clear error message if imports don't resolve. The README lists setup steps.

**Frame rate asymmetry in metrics** → Both systems are evaluated fairly at their native rates, but this means the raw TP/FP/FN counts aren't directly comparable across systems (ViT has 2x more frames). The F2/precision/recall *ratios* are still comparable. The README will note this.

**MPS compatibility** → Some PyTorch operations may not yet be supported on MPS. The script falls back to CPU automatically if MPS fails.

**Model checkpoint availability** → If HuggingFace is unreachable or the model path changes, the script fails at download. The checkpoint path is a constant that can be updated.
