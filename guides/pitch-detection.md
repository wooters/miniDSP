# Pitch Detection {#pitch-detection}

This [pitch detection](https://en.wikipedia.org/wiki/Pitch_detection_algorithm) tutorial compares two classic [fundamental-frequency (F0)](https://en.wikipedia.org/wiki/Fundamental_frequency) estimators:

- **[Autocorrelation](https://en.wikipedia.org/wiki/Autocorrelation) peak picking** (time domain)
- **[FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform) peak picking** (frequency domain)

Both are implemented in miniDSP and demonstrated in
`examples/pitch_detection.c`.

Build and run the example from the repository root:

```sh
make -C examples pitch_detection
cd examples && ./pitch_detection
open pitch_detection.html
```

---

## Autocorrelation F0

For a voiced frame, the fundamental period shows up as a strong peak in
the autocorrelation function:

\f[
R[\tau] = \frac{\sum_{n=0}^{N-1-\tau} x[n]x[n+\tau]}
               {\sum_{n=0}^{N-1} x[n]^2},
\qquad
f_0 = \frac{f_s}{\tau_{\text{peak}}}
\f]

We search only lags mapped from a desired F0 range
(`min_freq_hz..max_freq_hz`), then choose the strongest local peak.

### Reading the algorithm in C

```c
// x[n] -> frame[n], fs -> sample_rate
// lag_min/lag_max come from requested min/max F0 range.

double r0 = 0.0;  // denominator: SUM x[n]^2
for (unsigned n = 0; n < N; n++) {
    r0 += frame[n] * frame[n];
}

double best_r = -1.0;
unsigned best_lag = 0;

for (unsigned tau = lag_min; tau <= lag_max; tau++) {
    // numerator: SUM x[n] * x[n+tau]
    double num = 0.0;
    for (unsigned n = 0; n < N - tau; n++) {
        num += frame[n] * frame[n + tau];
    }

    double r_tau = (r0 > 0.0) ? (num / r0) : 0.0;  // normalised R[tau]

    // local-max check using neighbors (R[tau-1], R[tau], R[tau+1])
    if (r_tau > best_r /* and is_local_peak */) {
        best_r = r_tau;
        best_lag = tau;
    }
}

double f0_hz = (best_lag > 0) ? (sample_rate / (double)best_lag) : 0.0;
```

---

## FFT-based F0

This method applies a [Hann window](https://en.wikipedia.org/wiki/Hann_function), computes the one-sided FFT magnitude,
and picks the dominant peak in a frequency range:

\f[
f_0 = \frac{k_{\text{peak}} f_s}{N}
\f]

It is simple and fast, but more sensitive to noise and harmonic dominance
than autocorrelation.

### Reading the algorithm in C

```c
// x[n] -> frame[n], fs -> sample_rate
// 1) Apply Hann window:
for (unsigned n = 0; n < N; n++) {
    double w = 0.5 * (1.0 - cos(2.0 * M_PI * (double)n / (double)(N - 1)));
    xw[n] = frame[n] * w;
}

// 2) Compute magnitude spectrum directly from DFT definition (educational form):
for (unsigned k = 0; k <= N / 2; k++) {
    double re = 0.0, im = 0.0;
    for (unsigned n = 0; n < N; n++) {
        double phase = 2.0 * M_PI * (double)k * (double)n / (double)N;
        re += xw[n] * cos(phase);
        im -= xw[n] * sin(phase);
    }
    mag[k] = sqrt(re * re + im * im);
}

// 3) Search bins mapped from requested F0 range:
unsigned k_min = (unsigned)ceil(min_freq_hz * (double)N / sample_rate);
unsigned k_max = (unsigned)floor(max_freq_hz * (double)N / sample_rate);
unsigned k_peak = k_min;
for (unsigned k = k_min; k <= k_max; k++) {
    if (mag[k] > mag[k_peak]) k_peak = k;
}

double f0_hz = (double)k_peak * sample_rate / (double)N;
```

---

## Frame-Wise Tracking

In practice, pitch is estimated frame-by-frame over time:

\snippet pitch_detection.c frame-tracking

---

## Visual Comparison

\htmlonly
<p><strong>Ground truth vs estimated tracks (entire signal)</strong></p>
<iframe src="pitch_f0_tracks.html"
        style="width:100%;height:390px;border:1px solid #ddd;border-radius:4px;"
        frameborder="0"></iframe>

<p style="margin-top:0.8rem;"><strong>Autocorrelation peak (single frame)</strong></p>
<iframe src="pitch_acf_peak_frame.html"
        style="width:100%;height:390px;border:1px solid #ddd;border-radius:4px;"
        frameborder="0"></iframe>

<p style="margin-top:0.8rem;"><strong>FFT peak pick (single frame)</strong></p>
<iframe src="pitch_fft_peak_frame.html"
        style="width:100%;height:390px;border:1px solid #ddd;border-radius:4px;"
        frameborder="0"></iframe>
\endhtmlonly

---

## Failure Modes and Trade-offs

- **Autocorrelation** can fail on weakly voiced/noisy frames where no clear
  lag peak exists.
- **FFT peak pick** can lock onto [harmonics](https://en.wikipedia.org/wiki/Harmonic) (e.g. 2f0, 3f0) when the
  fundamental is weak.
- Restricting the search range (`min_freq_hz`, `max_freq_hz`) is critical for
  both methods.
- Short frames improve time resolution but reduce frequency/lag resolution.

Both miniDSP APIs return `0.0` when no reliable F0 peak is found.

---

## API Reference

- MD_f0_autocorrelation() -- F0 from autocorrelation peak
- MD_f0_fft() -- F0 from Hann-windowed FFT magnitude peak
- MD_shutdown() -- release cached FFT plans/resources when done
