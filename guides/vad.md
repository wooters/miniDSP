# Voice Activity Detection {#vad}

This tutorial builds a [voice activity detector](https://en.wikipedia.org/wiki/Voice_activity_detection) (VAD) that combines five normalized audio features into a weighted score, applies a threshold, and smooths the binary decision with onset gating and hangover.  VAD is a fundamental preprocessing step for streaming audio, [ASR](https://en.wikipedia.org/wiki/Speech_recognition) pipelines, and noise reduction.

All features are computed from existing miniDSP primitives.  The detector adapts to different microphones, gain settings, and noise floors without manual tuning.

Build and run the example from the repository root:

```sh
make -C examples vad
cd examples && ./vad
open vad_plot.html
```

---

## Energy

Frame energy measures overall loudness — the simplest voice activity cue.  Speech frames carry more energy than silence or low-level background noise.

\f[
E = \sum_{n=0}^{N-1} x[n]^2
\f]

**Reading the formula in C:**

```c
// x[n] -> signal[n], N -> frame length
double energy = 0.0;
for (unsigned n = 0; n < N; n++) {
    energy += signal[n] * signal[n];  // x[n]^2
}
```

**Intuition:** silence has near-zero energy; speech has high energy.  Energy alone fails in moderate noise because the noise floor raises the baseline.

---

## Zero-crossing rate

The [zero-crossing rate](https://en.wikipedia.org/wiki/Zero-crossing_rate) counts how often the signal changes sign, normalized by the number of adjacent pairs:

\f[
\mathrm{ZCR} = \frac{1}{N-1}\sum_{n=1}^{N-1}
    \mathbf{1}\!\bigl[\mathrm{sgn}(x[n]) \ne \mathrm{sgn}(x[n-1])\bigr]
\f]

**Reading the formula in C:**

```c
// x[n] -> signal[n], N -> frame length
unsigned crossings = 0;
for (unsigned n = 1; n < N; n++) {
    // 1[sgn(x[n]) != sgn(x[n-1])]
    if ((signal[n] >= 0.0) != (signal[n - 1] >= 0.0))
        crossings++;
}
double zcr = (double)crossings / (double)(N - 1);  // normalize to [0, 1]
```

**Intuition:** voiced speech has a low ZCR (periodic, low-frequency fundamental); unvoiced fricatives have high ZCR; silence has low ZCR.  ZCR helps distinguish voiced speech from broadband noise.

---

## Spectral entropy

[Spectral entropy](https://en.wikipedia.org/wiki/Spectral_density#Spectral_entropy) measures how uniformly energy is spread across frequency bins.  Low entropy means energy is concentrated (tonal); high entropy means energy is diffuse (noise-like).

\f[
H = -\frac{1}{\ln(K)} \sum_{k=0}^{K-1} p_k \ln(p_k),
\qquad p_k = \frac{\mathrm{PSD}[k]}{\sum_{j=0}^{K-1} \mathrm{PSD}[j]}
\f]

where \f$K = N/2 + 1\f$ is the number of one-sided PSD bins.

**Reading the formula in C:**

```c
// PSD[k] -> psd[k], K -> num_bins
double total = 0.0;
for (unsigned k = 0; k < num_bins; k++)
    total += psd[k];  // denominator of p_k

double entropy = 0.0;
for (unsigned k = 0; k < num_bins; k++) {
    double p_k = psd[k] / total;         // p_k = PSD[k] / sum(PSD)
    if (p_k > 0.0)
        entropy -= p_k * log(p_k);       // -sum(p_k * ln(p_k))
}
entropy /= log((double)num_bins);         // normalize by ln(K) -> [0, 1]
```

**Intuition:** speech has lower spectral entropy (energy concentrated at harmonics); noise has higher spectral entropy (energy spread uniformly).

---

## Spectral flatness

[Spectral flatness](https://en.wikipedia.org/wiki/Spectral_flatness) (also called the Wiener entropy) is the ratio of the geometric mean to the arithmetic mean of PSD bins:

\f[
\mathrm{SF} = \frac{\left(\prod_{k=0}^{K-1} \mathrm{PSD}[k]\right)^{1/K}}
                    {\frac{1}{K}\sum_{k=0}^{K-1} \mathrm{PSD}[k]}
\f]

**Reading the formula in C:**

```c
// PSD[k] -> psd[k], K -> num_bins
// Compute in log domain to avoid overflow in the product
double log_sum = 0.0;
double arith_sum = 0.0;
for (unsigned k = 0; k < num_bins; k++) {
    double val = psd[k] > 0.0 ? psd[k] : 1e-30;
    log_sum += log(val);         // sum of ln(PSD[k])
    arith_sum += psd[k];        // sum of PSD[k]
}
double log_geo_mean = log_sum / (double)num_bins;  // (1/K) * sum(ln(PSD[k]))
double geo_mean = exp(log_geo_mean);               // geometric mean
double arith_mean = arith_sum / (double)num_bins;   // arithmetic mean

double flatness = geo_mean / arith_mean;  // SF in [0, 1]
```

**Intuition:** SF = 1.0 for white noise (perfectly flat spectrum); SF approaches 0 for a pure tone (all energy in one bin).  Speech falls between these extremes, with lower flatness than background noise.

---

## Band energy ratio

The band energy ratio measures the concentration of energy in the speech band (typically 300–3400 Hz, the [telephone bandwidth](https://en.wikipedia.org/wiki/Voice_frequency)):

\f[
\mathrm{BER} = \frac{\sum_{k : f_k \in [f_{\mathrm{lo}},\, f_{\mathrm{hi}}]} \mathrm{PSD}[k]}
                     {\sum_{k=0}^{K-1} \mathrm{PSD}[k]}
\f]

where \f$f_k = k \cdot f_s / N\f$ is the frequency of bin \f$k\f$.

**Reading the formula in C:**

```c
// PSD[k] -> psd[k], K -> num_bins
// f_lo -> band_low_hz, f_hi -> band_high_hz
double freq_per_bin = sample_rate / (double)N;  // f_s / N
double total = 0.0;
double band = 0.0;
for (unsigned k = 0; k < num_bins; k++) {
    double freq = k * freq_per_bin;  // f_k
    total += psd[k];
    if (freq >= band_low_hz && freq <= band_high_hz)
        band += psd[k];  // sum PSD in speech band
}
double ber = band / total;  // BER in [0, 1]
```

**Intuition:** speech concentrates energy in the 300–3400 Hz band.  Background noise typically has a flatter distribution, yielding a lower BER.

---

## Adaptive normalization

Raw feature values vary by orders of magnitude across microphones and gain settings.  The VAD normalizes each feature to [0, 1] using per-feature min/max estimates tracked via [exponential moving average](https://en.wikipedia.org/wiki/Moving_average#Exponential_moving_average) (EMA):

\f[
m_i \leftarrow m_i + \alpha \cdot (f_i - m_i), \quad
M_i \leftarrow M_i + \alpha \cdot (f_i - M_i)
\f]

where \f$m_i\f$ and \f$M_i\f$ are the running minimum and maximum for feature \f$i\f$, \f$f_i\f$ is the current raw value, and \f$\alpha\f$ is the adaptation rate.  The normalized feature is:

\f[
\hat{f}_i = \frac{f_i - m_i}{M_i - m_i}
\f]

clamped to [0, 1].

**Reading the formula in C:**

```c
// alpha -> adaptation_rate, f_i -> raw[i]
// m_i -> feat_min[i], M_i -> feat_max[i]
// Update min/max via EMA (only when raw value exceeds current bound)
if (raw[i] < feat_min[i])
    feat_min[i] = feat_min[i] + alpha * (raw[i] - feat_min[i]);
if (raw[i] > feat_max[i])
    feat_max[i] = feat_max[i] + alpha * (raw[i] - feat_max[i]);

// Normalize to [0, 1]
double range = feat_max[i] - feat_min[i];
if (range < 1e-12) range = 1e-12;  // prevent division by zero
double norm = (raw[i] - feat_min[i]) / range;
if (norm < 0.0) norm = 0.0;
if (norm > 1.0) norm = 1.0;
```

---

## Weighted scoring

The five normalized features are combined into a single score via weighted sum:

\f[
S = \sum_{i=0}^{4} w_i \cdot \hat{f}_i
\f]

where \f$w_i\f$ are caller-tunable weights (default: 0.2 each for equal weighting).

**Reading the formula in C:**

```c
// w_i -> weights[i], f_hat_i -> norm[i]
double score = 0.0;
for (int i = 0; i < 5; i++) {
    score += weights[i] * norm[i];  // S = sum(w_i * f_hat_i)
}
```

**Intuition:** equal weights give each feature equal vote.  Setting one weight to 1.0 and the rest to 0.0 reduces the detector to a single-feature VAD.

---

## State machine

The detector uses an onset + hangover mechanism to smooth the binary decision:

| State       | Condition                              | Action                           |
|:------------|:---------------------------------------|:---------------------------------|
| **SILENCE** | score >= threshold                     | Increment onset counter          |
| **SILENCE** | onset counter >= onset_frames          | Transition to SPEECH             |
| **SILENCE** | score < threshold                      | Reset onset counter              |
| **SPEECH**  | score >= threshold                     | Reset hangover to hangover_frames|
| **SPEECH**  | score < threshold                      | Decrement hangover counter       |
| **SPEECH**  | hangover counter == 0                  | Transition to SILENCE            |

**Onset gating** prevents transient clicks from triggering false positives — the score must exceed the threshold for `onset_frames` consecutive frames before the detector declares speech.

**Hangover** bridges brief dips mid-utterance — after the score drops below threshold, the speech decision holds for `hangover_frames` additional frames.

---

## Visualization

The interactive plot below shows the VAD processing a synthesized signal with two speech segments (sine bursts at 1000 Hz) separated by silence:

\htmlonly
<iframe src="vad_plot.html" width="100%" height="850" frameborder="0"
        style="border: 1px solid #e2e8f0; border-radius: 6px;"></iframe>
\endhtmlonly

The four panels show:
1. **Waveform** — peak envelope per frame, showing speech and silence regions.
2. **Normalized features** — all five features mapped to [0, 1], showing how they respond differently to speech vs. silence.
3. **Combined score** — the weighted sum of features compared against the threshold (dashed red line).
4. **Decision** — the final binary output after onset gating and hangover smoothing.

---

## API summary

**API:**

- `MD_vad_default_params(MD_vad_params *params)` — populate with sensible defaults.
- `MD_vad_init(MD_vad_state *state, const MD_vad_params *params)` — initialize state (NULL params for defaults).
- `MD_vad_calibrate(MD_vad_state *state, const double *signal, unsigned N, double sample_rate)` — feed known-silence frames to seed normalization.
- `MD_vad_process_frame(MD_vad_state *state, const double *signal, unsigned N, double sample_rate, double *score_out, double *features_out)` — process one frame, return 0 (silence) or 1 (speech).

**Quick example** — initialize and process:

\snippet vad.c vad-init

\snippet vad.c vad-process

**Calibration** — optional, improves accuracy if silence frames are available:

\snippet vad.c vad-calibrate

**Custom weights** — use a single feature:

\snippet vad.c vad-custom-weights
