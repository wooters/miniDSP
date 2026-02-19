# Mel Filterbank and MFCCs {#mel-mfcc}

This tutorial introduces two classic speech/audio front-end features:

- **Mel filterbank energies** — triangular spectral bands spaced on the
  [mel scale](https://en.wikipedia.org/wiki/Mel_scale).
- **[MFCCs](https://en.wikipedia.org/wiki/Mel-frequency_cepstrum)** — DCT of log mel energies.

miniDSP provides both as frame-level APIs in `minidsp.h`.

## Why mel and MFCC?

The linear FFT axis over-resolves high frequencies compared to human pitch
perception. Mel filterbanks compress frequency spacing to be denser at low
frequencies and coarser at high frequencies.

MFCCs then decorrelate log mel energies via a DCT, producing compact features
used widely in speech recognition and audio classification.

## Step 1: Build mel energies

miniDSP uses:

- HTK mel mapping:
  \f[
  m(f) = 2595 \log_{10}\!\left(1 + \frac{f}{700}\right)
  \f]
- Internal Hann windowing
- One-sided PSD bins: \f$|X(k)|^2 / N\f$

Compute one frame of mel energies:

\snippet mel_mfcc.c compute-mel

## Step 2: Compute MFCCs

MFCCs are computed from log mel energies with a DCT-II:

\f[
c_n = \alpha_n \sum_{m=0}^{M-1}
\log(\max(E_m, 10^{-12}))
\cos\!\left(\frac{\pi n (m + 1/2)}{M}\right)
\f]

where:

- \f$M\f$ is the number of mel bands
- \f$\alpha_0 = \sqrt{1/M}\f$
- \f$\alpha_n = \sqrt{2/M}\f$ for \f$n > 0\f$

miniDSP returns **C0 in `mfcc_out[0]`**.

\snippet mel_mfcc.c compute-mfcc

## Practical notes

- Requested frequency bounds are runtime-clamped to `[0, Nyquist]`.
- If the clamped band is empty, mel energies are zero and MFCCs remain finite
  via the log floor.
- `MD_shutdown()` should be called when done with FFT-based APIs.

## API reference

- MD_mel_filterbank() — build mel triangular weight matrix
- MD_mel_energies() — compute mel-band energies from one frame
- MD_mfcc() — compute MFCC vector (C0 included)
- MD_shutdown() — release cached FFT resources

