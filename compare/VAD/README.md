# VAD system comparison

Compare miniDSP VAD against the ViT-MFCC baseline from the
[LibriVAD](https://arxiv.org/abs/2512.17281) project on the same test data.

## Systems

| System | Type | Parameters | Frame rate |
|--------|------|-----------|------------|
| miniDSP VAD | Hand-crafted features + weighted threshold + state machine | ~10 tunable hyperparams | 20 ms |
| ViT-MFCC (small) | Vision Transformer on 39-dim MFCCs | ~3.5M trained weights | 10 ms |

The ViT-MFCC model is the KWT (Keyword Transformer) architecture trained on
LibriSpeechConcat with cross-entropy loss. The pre-trained checkpoint is
downloaded automatically from HuggingFace.

## Prerequisites

- [uv](https://docs.astral.sh/uv/getting-started/installation/)
- The [LibriVAD](https://github.com/LibriVAD/LibriVAD) repo with test data
  prepared (noisy WAVs in `Results/` and labels in `Files/Labels/`)

```bash
export LIBRIVAD_DIR=~/projects/LibriVAD
```

## Setup

```bash
cd compare/VAD
uv sync
```

## Usage

Full comparison on test-clean with per-condition breakdown:

```bash
uv run python compare_vad.py \
    --librivad-root ${LIBRIVAD_DIR} \
    --breakdown
```

Quick sanity check on a subset:

```bash
uv run python compare_vad.py \
    --librivad-root ${LIBRIVAD_DIR} \
    --noises Babble_noise SSN_noise \
    --snrs 10 20
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `--librivad-root` | *(required)* | Path to LibriVAD project root |
| `--dataset` | `LibriSpeechConcat` | Dataset variant |
| `--split` | `test-clean` | Dataset split |
| `--noises` | all | Noise types to include |
| `--snrs` | all | SNR levels to include |
| `--beta` | `2.0` | F-beta parameter (>1 favors recall) |
| `--breakdown` | off | Print per-noise and per-SNR tables |

## Interpreting results

### Frame rate asymmetry

The two systems operate at different frame rates: miniDSP at 20 ms and
ViT-MFCC at 10 ms. Each system is evaluated against ground-truth labels
downsampled to its own frame rate via majority voting. This means:

- The ViT has 2x finer temporal resolution, which is an inherent design
  advantage
- Raw TP/FP/FN counts are not directly comparable across systems
- The F2/precision/recall *ratios* are comparable and are the primary metrics

### Metric definitions

- **Precision**: TP / (TP + FP) -- of frames labeled speech, how many were
  actually speech
- **Recall**: TP / (TP + FN) -- of actual speech frames, how many were
  detected
- **F2**: Weighted harmonic mean with beta=2, giving recall twice the weight
  of precision. F2 = 5PR / (4P + R)

### ViT threshold

The ViT uses a fixed 0.5 threshold (argmax of 2-class softmax), which is the
designed operating point for a cross-entropy-trained classifier. The miniDSP
threshold was optimized via Optuna (0.245), giving it a tuned operating point.
