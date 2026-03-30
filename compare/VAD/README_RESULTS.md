# VAD Comparison Results

Comparison of miniDSP VAD against the ViT-MFCC (small) baseline on the
LibriVAD **test-clean** split (LibriSpeechConcat dataset).

**Date:** 2026-03-30
**Hardware:** Mac Mini, Apple Silicon (MPS backend for PyTorch)

## Systems Under Test

| Property | miniDSP VAD | ViT-MFCC (small) |
|----------|-------------|-------------------|
| Type | Hand-crafted features + weighted threshold + onset/hangover state machine | Vision Transformer (KWT) on MFCCs |
| Language | C (via pyminidsp CFFI bindings) | Python / PyTorch |
| Trainable parameters | 0 (10 tuned hyperparameters) | ~3.5M |
| Frame rate | 20 ms | 10 ms |
| Training data | Optuna optimization on train-clean-100 | Supervised training on LibriSpeechConcat (small) |
| Threshold | 0.245 (optimized) | 0.5 (argmax, designed operating point) |
| Features | Energy, ZCR, spectral entropy, spectral flatness, band energy ratio | 39-dim MFCCs (12 + energy + delta + delta-delta) |

## Overall Results

702 files (9 noise types x 6 SNR levels x ~13 files each; 54 files skipped due to missing labels).

| System | F2 | Precision | Recall | Wall time |
|--------|---:|----------:|-------:|----------:|
| miniDSP VAD | 0.8440 | 0.8274 | 0.8482 | 4.3 s |
| ViT-MFCC (small) | 0.9614 | 0.9390 | 0.9672 | 79.3 s |

The ViT-MFCC leads by **+0.117 F2** but takes **18x longer** to process the
same data.

## Per SNR Breakdown

| SNR (dB) | miniDSP F2 | ViT F2 | Gap |
|---------:|-----------:|-------:|----:|
| -5 | 0.7993 | 0.9321 | 0.133 |
| 0 | 0.8358 | 0.9523 | 0.117 |
| 5 | 0.8512 | 0.9606 | 0.109 |
| 10 | 0.8515 | 0.9685 | 0.117 |
| 15 | 0.8588 | 0.9752 | 0.116 |
| 20 | 0.8688 | 0.9801 | 0.111 |

The gap is relatively consistent across SNR levels (~0.11-0.13), widening
slightly at the hardest condition (-5 dB).

## Per Noise Type Breakdown

| Noise Type | miniDSP F2 | ViT F2 | Gap |
|------------|------------|--------|-----|
| Babble_noise | 0.8792 | 0.9289 | 0.050 |
| City_noise | 0.8553 | 0.9523 | 0.097 |
| Domestic_noise | 0.8524 | 0.9710 | 0.119 |
| Nature_noise | 0.8252 | 0.9761 | 0.151 |
| Office_noise | 0.8069 | 0.9760 | 0.169 |
| Public_noise | 0.8595 | 0.9485 | 0.089 |
| SSN_noise | 0.8827 | 0.9623 | 0.080 |
| Street_noise | 0.8291 | 0.9618 | 0.133 |
| Transport_noise | 0.7996 | 0.9773 | 0.178 |

miniDSP does best on noise types with distinctive spectral characteristics
(babble, SSN) where its energy and band-energy-ratio features are most
discriminative. The largest gaps appear with transport and office noise, which
overlap more with the speech band.

## Key Takeaways

1. **miniDSP VAD achieves 88% of the ViT's F2 score** (0.844 / 0.961) with
   zero trainable parameters, pure C implementation, and 18x faster inference.

2. **The speed difference is significant**: 4.3 s vs 79.3 s for 702 files.
   The miniDSP VAD is designed for real-time embedded use; the ViT requires
   a GPU (or Apple MPS) for practical throughput.

3. **The miniDSP recall drops on test-clean** compared to train-clean-100
   (0.848 vs 0.981), suggesting the Optuna optimization overfit somewhat to
   the training conditions. The ViT generalizes better from its training data.

4. **Noise-type sensitivity differs**: miniDSP struggles most with transport
   and office noise (spectral overlap with speech), while the ViT is more
   robust across all noise types.

5. **For applications where latency, footprint, and simplicity matter**
   (embedded devices, real-time pipelines), the miniDSP VAD offers a strong
   cost/performance tradeoff. For maximum accuracy with compute budget
   available, the ViT-MFCC is clearly superior.

## Methodology Notes

- Each system is evaluated at its native frame rate (miniDSP at 20 ms, ViT at
  10 ms). Ground-truth labels are downsampled independently to match.
- F2 metric (beta=2) weights recall twice as heavily as precision.
- The ViT uses a fixed 0.5 decision threshold (argmax of 2-class softmax).
- The miniDSP VAD uses its Optuna-optimized defaults from the C library.
- Wall times include audio I/O and feature extraction, not model loading.
