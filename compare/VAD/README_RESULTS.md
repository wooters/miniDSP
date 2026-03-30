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

| System | F2 | Precision | Recall | AUC (macro) | AUC (pooled) | Wall time |
|--------|---:|----------:|-------:|------------:|-------------:|----------:|
| miniDSP VAD | 0.8440 | 0.8274 | 0.8482 | 0.6519 | 0.6517 | 3.6 s |
| ViT-MFCC (small) | 0.9614 | 0.9390 | 0.9672 | 0.9712 | 0.9819 | 77.8 s |

AUC (macro) = mean of per-file AUC values; AUC (pooled) = AUC on
concatenated frames. Macro-averaged AUC matches the methodology used in the
[LibriVAD paper](https://arxiv.org/abs/2512.17281) (Table 7), which reports
0.9710 for MFCC-ViT on this dataset -- our 0.9712 confirms reproducibility.

The ViT-MFCC leads by **+0.117 F2** and **+0.319 AUC (macro)** but takes
**22x longer** to process the same data.

## Per SNR Breakdown

| SNR (dB) | miniDSP F2 | ViT F2 | F2 Gap | miniDSP AUC | ViT AUC | AUC Gap |
|---------:|-----------:|-------:|-------:|------------:|--------:|--------:|
| -5 | 0.7993 | 0.9321 | 0.133 | 0.5785 | 0.9145 | 0.336 |
| 0 | 0.8358 | 0.9523 | 0.117 | 0.6205 | 0.9564 | 0.336 |
| 5 | 0.8512 | 0.9606 | 0.109 | 0.6539 | 0.9787 | 0.325 |
| 10 | 0.8515 | 0.9685 | 0.117 | 0.6730 | 0.9888 | 0.316 |
| 15 | 0.8588 | 0.9752 | 0.116 | 0.6870 | 0.9934 | 0.306 |
| 20 | 0.8688 | 0.9801 | 0.111 | 0.6985 | 0.9955 | 0.297 |

AUC values are macro-averaged per file within each SNR bucket.

The F2 gap is relatively consistent across SNR levels (~0.11-0.13), widening
slightly at the hardest condition (-5 dB). The AUC gap is larger (~0.30-0.34)
and narrows with increasing SNR, showing that miniDSP's continuous scores
improve more with cleaner signals than its binary decisions suggest.

## Per Noise Type Breakdown

| Noise Type | miniDSP F2 | ViT F2 | F2 Gap | miniDSP AUC | ViT AUC | AUC Gap |
|------------|------------|--------|-------:|------------:|--------:|--------:|
| Babble_noise | 0.8792 | 0.9289 | 0.050 | 0.5970 | 0.8556 | 0.259 |
| City_noise | 0.8553 | 0.9523 | 0.097 | 0.6499 | 0.9760 | 0.326 |
| Domestic_noise | 0.8524 | 0.9710 | 0.119 | 0.7172 | 0.9912 | 0.274 |
| Nature_noise | 0.8252 | 0.9761 | 0.151 | 0.6644 | 0.9955 | 0.331 |
| Office_noise | 0.8069 | 0.9760 | 0.169 | 0.6751 | 0.9957 | 0.321 |
| Public_noise | 0.8595 | 0.9485 | 0.089 | 0.6356 | 0.9642 | 0.329 |
| SSN_noise | 0.8827 | 0.9623 | 0.080 | 0.6180 | 0.9858 | 0.368 |
| Street_noise | 0.8291 | 0.9618 | 0.133 | 0.6467 | 0.9823 | 0.336 |
| Transport_noise | 0.7996 | 0.9773 | 0.178 | 0.6634 | 0.9946 | 0.331 |

AUC values are macro-averaged per file within each noise bucket.

miniDSP does best on noise types with distinctive spectral characteristics
(babble, SSN) where its energy and band-energy-ratio features are most
discriminative. The largest F2 gaps appear with transport and office noise,
which overlap more with the speech band. AUC tells a different story:
miniDSP's best AUC is on domestic noise (0.717), suggesting its continuous
scores are most calibrated for that category, while babble and SSN -- where
F2 is strongest -- have relatively low AUC (0.60-0.62), indicating the good
F2 relies heavily on threshold tuning rather than score quality.

## Key Takeaways

1. **miniDSP VAD achieves 88% of the ViT's F2 score** (0.844 / 0.961) with
   zero trainable parameters, pure C implementation, and 19x faster inference.

2. **AUC reveals a larger quality gap than F2**: miniDSP macro-averaged AUC
   is 0.652 vs ViT's 0.971 -- a +0.319 gap compared to +0.117 for F2. This
   means the miniDSP's continuous scores have limited discriminative power;
   its good F2 relies on threshold tuning by Optuna rather than intrinsic
   score quality. Our ViT AUC of 0.9712 closely reproduces the 0.9710
   reported in the [LibriVAD paper](https://arxiv.org/abs/2512.17281)
   (Table 7).

3. **The speed difference is significant**: 3.6 s vs 77.8 s for 702 files.
   The miniDSP VAD is designed for real-time embedded use; the ViT requires
   a GPU (or Apple MPS) for practical throughput.

4. **The miniDSP recall drops on test-clean** compared to train-clean-100
   (0.848 vs 0.981), suggesting the Optuna optimization overfit somewhat to
   the training conditions. The ViT generalizes better from its training data.

5. **Noise-type sensitivity differs**: miniDSP struggles most with transport
   and office noise (spectral overlap with speech), while the ViT is more
   robust across all noise types.

6. **For applications where latency, footprint, and simplicity matter**
   (embedded devices, real-time pipelines), the miniDSP VAD offers a strong
   cost/performance tradeoff. For maximum accuracy with compute budget
   available, the ViT-MFCC is clearly superior. The low AUC suggests that
   if the miniDSP threshold needs to be re-tuned for different conditions,
   there is limited headroom in the score distribution.

## Methodology Notes

- Each system is evaluated at its native frame rate (miniDSP at 20 ms, ViT at
  10 ms). Ground-truth labels are downsampled independently to match.
- F2 metric (beta=2) weights recall twice as heavily as precision.
- The ViT uses a fixed 0.5 decision threshold (argmax of 2-class softmax).
- The miniDSP VAD uses its Optuna-optimized defaults from the C library.
- Wall times include audio I/O and feature extraction, not model loading.
- **AUC (macro)**: AUC-ROC is computed per file, then averaged across files.
  This matches the LibriVAD paper's methodology and gives equal weight to
  every file regardless of length.
- **AUC (pooled)**: AUC-ROC is computed on all frames concatenated. Longer
  files contribute more frames and therefore more weight. Pooled AUC tends to
  be higher because easy/long files dominate.
- miniDSP continuous scores are the weighted feature combination (0-1) from
  the C library. ViT scores are softmax speech probabilities.
