## ADDED Requirements

### Requirement: Compare script evaluates both VAD systems on LibriVAD test data
The script SHALL load audio files from a LibriVAD Results directory, run both miniDSP VAD and ViT-MFCC inference on each file, and compute F2, precision, and recall for both systems against ground-truth labels.

#### Scenario: Full test-clean evaluation
- **WHEN** the user runs `compare_vad.py` with `--librivad-root` pointing to a LibriVAD project containing prepared test-clean data
- **THEN** the script evaluates all 756 files (9 noise types × 6 SNRs × 14 files) and prints overall F2, precision, and recall for both miniDSP VAD and ViT-MFCC

#### Scenario: Filtered evaluation
- **WHEN** the user specifies `--noises` and/or `--snrs` to restrict conditions
- **THEN** only matching files are evaluated and metrics are computed on the subset

### Requirement: Script auto-downloads ViT-MFCC checkpoint from HuggingFace
The script SHALL download the pre-trained ViT-MFCC (small, LibriSpeechConcat) model checkpoint from HuggingFace using `huggingface_hub.hf_hub_download` if not already cached locally.

#### Scenario: First run downloads model
- **WHEN** the script is run for the first time and no cached checkpoint exists
- **THEN** the checkpoint is downloaded from `LibriVAD/LibriVAD` dataset repo at path `Eval-ViT-MFCC/small-models/LibriSpeechConcat_vit_MFCC/model_epoch_50.model` and cached for future runs

#### Scenario: Subsequent runs use cache
- **WHEN** the script is run again after a prior download
- **THEN** the cached checkpoint is loaded without re-downloading

### Requirement: miniDSP VAD evaluated with library defaults
The script SHALL instantiate `pyminidsp.VAD()` with no arguments (using the Optuna-optimized C library defaults) and call `.process()` at 20ms frame rate.

#### Scenario: Default parameters used
- **WHEN** the miniDSP VAD is evaluated
- **THEN** it uses the default parameters from the C library (threshold=0.245332, onset_frames=1, hangover_frames=22, etc.) without any overrides

### Requirement: ViT-MFCC evaluated at 0.5 threshold with reference code
The script SHALL import the KWT model class and MFCC extraction from the LibriVAD `Eval-ViT-MFCC/` directory, extract 39-dim MFCC features at 10ms frame rate, run inference through the pre-trained model, and apply argmax (threshold 0.5) for binary decisions.

#### Scenario: ViT inference pipeline
- **WHEN** a wav file is processed through the ViT-MFCC pipeline
- **THEN** MFCCs are extracted (12 coefficients + energy + delta + delta-delta = 39 dims, 25ms window, 10ms shift), CMVN-normalized, split into 100-frame sequences, passed through the KWT model, softmax applied, and frames with P(speech) >= 0.5 are classified as speech

### Requirement: Labels downsampled to each system's native frame rate
The script SHALL downsample sample-level ground-truth labels to frame-level via majority voting, independently for each system's frame rate (20ms for miniDSP, 10ms for ViT-MFCC).

#### Scenario: Label alignment
- **WHEN** labels are prepared for metric computation
- **THEN** miniDSP labels are downsampled at 320 samples/frame (20ms × 16kHz) and ViT labels are downsampled at 160 samples/frame (10ms × 16kHz), both using majority vote per frame

### Requirement: Per-condition breakdown output
The script SHALL support a `--breakdown` flag that prints metrics grouped by noise type and by SNR level for both systems.

#### Scenario: Breakdown by noise and SNR
- **WHEN** `--breakdown` is specified
- **THEN** the output includes a per-noise-type table and a per-SNR table showing F2, precision, and recall for both systems side by side

### Requirement: Apple Silicon GPU acceleration
The script SHALL detect and use Apple Silicon MPS acceleration for PyTorch inference when available, falling back to CUDA then CPU.

#### Scenario: MPS available
- **WHEN** running on Apple Silicon with MPS-capable PyTorch
- **THEN** ViT inference uses the MPS device

#### Scenario: MPS unavailable
- **WHEN** MPS is not available
- **THEN** the script falls back to CUDA if available, otherwise CPU

### Requirement: Self-contained uv project with README
The `compare/VAD/` directory SHALL contain a `pyproject.toml` for dependency management via `uv` and a `README.md` documenting setup, usage, prerequisites, and how to interpret results.

#### Scenario: Clean setup
- **WHEN** a user clones the repo and navigates to `compare/VAD/`
- **THEN** running `uv sync` installs all dependencies and the script is ready to run
