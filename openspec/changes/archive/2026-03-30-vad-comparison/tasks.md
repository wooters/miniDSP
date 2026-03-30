## 1. Project Setup

- [x] 1.1 Create `compare/VAD/` directory with `pyproject.toml` (uv project) listing all dependencies: torch, huggingface_hub, einops, scipy, numpy, scikit-learn, soundfile, pyminidsp
- [x] 1.2 Create `compare/VAD/README.md` documenting prerequisites (LibriVAD repo with test-clean data prepared), setup (`uv sync`), usage examples, and notes on frame-rate asymmetry and metric interpretation

## 2. Core Comparison Script

- [x] 2.1 Implement CLI argument parsing: `--librivad-root`, `--dataset` (default LibriSpeechConcat), `--split` (default test-clean), `--noises`, `--snrs`, `--beta` (default 2.0), `--breakdown`, `--workers`
- [x] 2.2 Implement LibriVAD file discovery and label loading (reuse logic from `optimize_vad.py`: `resolve_label_path`, `sample_labels_to_frame_labels`, `discover_librivad_files`, `load_librivad_dataset`)
- [x] 2.3 Implement ViT-MFCC model download via `huggingface_hub.hf_hub_download` (repo_id=`LibriVAD/LibriVAD`, repo_type=`dataset`, filename=`Eval-ViT-MFCC/small-models/LibriSpeechConcat_vit_MFCC/model_epoch_50.model`)
- [x] 2.4 Implement ViT-MFCC inference pipeline: add `{librivad_root}/Eval-ViT-MFCC` to sys.path, import KWT/mfcc/cmvn/split_data_into_sequences/enframe, extract 39-dim MFCCs at 10ms frame rate, CMVN normalize, chunk into 100-frame sequences, forward pass with MPS/CUDA/CPU device selection, softmax, threshold at 0.5
- [x] 2.5 Implement miniDSP VAD evaluation: instantiate `VAD()` with defaults, call `.process()` at 20ms frame rate, collect binary decisions
- [x] 2.6 Implement F-beta metric computation (copy `compute_fbeta` logic from `optimize_vad.py`): accumulate TP/FP/FN per file, compute precision/recall/F2 overall and per condition

## 3. Output and Reporting

- [x] 3.1 Implement overall results table: side-by-side F2, precision, recall for both systems
- [x] 3.2 Implement `--breakdown` output: per-noise-type table and per-SNR table with both systems' metrics side by side

## 4. Validation

- [x] 4.1 Run the comparison script on test-clean and verify both systems produce plausible metrics (miniDSP F2 in the ~0.93 range based on train-clean-100 results; ViT produces non-trivial predictions)
- [x] 4.2 Verify HuggingFace download works and model loads correctly on Apple Silicon (MPS or CPU fallback)
