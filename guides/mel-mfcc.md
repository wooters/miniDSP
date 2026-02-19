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

**Reading the formula in C:**

```c
// f -> freq_hz, m(f) -> mel
double mel = 2595.0 * log10(1.0 + freq_hz / 700.0);
```

Compute one frame of mel energies:

\snippet mel_mfcc.c compute-mel

## Step 2: Compute MFCCs

MFCCs are computed from log mel energies with a DCT-II:

\f[
c_n = \alpha_n \sum_{m=0}^{M-1}
\log(\max(E_m, 10^{-12}))
\cos\!\left(\frac{\pi n (m + 1/2)}{M}\right)
\f]

**Reading the formula in C:**

```c
// E_m -> mel_energy[m], c_n -> mfcc[n], M -> num_mels, n -> coeff index
for (unsigned n = 0; n < num_coeffs; n++) {
    double alpha = (n == 0)
        ? sqrt(1.0 / (double)num_mels)
        : sqrt(2.0 / (double)num_mels);

    double acc = 0.0;
    for (unsigned m = 0; m < num_mels; m++) {
        double log_em = log(fmax(mel_energy[m], 1e-12));
        double basis = cos(M_PI * (double)n * ((double)m + 0.5) / (double)num_mels);
        acc += log_em * basis;
    }
    mfcc[n] = alpha * acc;
}
```

where:

- \f$M\f$ is the number of mel bands
- \f$\alpha_0 = \sqrt{1/M}\f$
- \f$\alpha_n = \sqrt{2/M}\f$ for \f$n > 0\f$

miniDSP returns **C0 in `mfcc_out[0]`**.

\snippet mel_mfcc.c compute-mfcc

## Visualisations

These plots are generated from one deterministic signal:
\f$x[n] = 0.7\sin(2\pi\cdot440t) + 0.2\cos(2\pi\cdot1000t) + 0.1\sin(2\pi\cdot3000t)\f$,
with \f$t=n/f_s\f$ and \f$f_s=8192\f$ Hz.
Mel energies and MFCCs are computed from the first 1024-sample analysis frame.

\htmlonly
<div style="display:flex;gap:12px;flex-wrap:wrap;">
  <iframe src="mel_input_waveform.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="mel_input_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
<div style="display:flex;gap:12px;flex-wrap:wrap;margin-top:0.8rem;">
  <iframe src="mel_filterbank_shapes.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="mel_energies_frame.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
<p style="margin-top:0.8rem;">
  <iframe src="mfcc_frame.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</p>
\endhtmlonly

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
