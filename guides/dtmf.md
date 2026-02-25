# DTMF Tone Detection and Generation {#dtmf}

[Dual-Tone Multi-Frequency (DTMF)](https://en.wikipedia.org/wiki/Dual-tone_multi-frequency_signaling)
is the signalling system used by touch-tone telephones.
Each keypad button is encoded as the sum of two sinusoids -- one from
a low-frequency **row** group and one from a high-frequency **column** group.
The receiver decodes the button by identifying both frequencies.

miniDSP provides ITU-T Q.24-compliant detection and generation in
`src/minidsp_dtmf.c`, demonstrated in `examples/dtmf_detector.c`.

Build and run the self-test from the repository root:

```sh
make -C examples dtmf_detector
cd examples && ./dtmf_detector
```

---

## The DTMF frequency table

Each button sits at the intersection of one row and one column frequency:

|            | 1209 Hz | 1336 Hz | 1477 Hz | 1633 Hz |
|:----------:|:-------:|:-------:|:-------:|:-------:|
| **697 Hz** |    1    |    2    |    3    |    A    |
| **770 Hz** |    4    |    5    |    6    |    B    |
| **852 Hz** |    7    |    8    |    9    |    C    |
| **941 Hz** |    *    |    0    |    #    |    D    |

The frequencies were chosen so that no tone is a harmonic of another
(ratios are never simple integers), preventing false triggers from
harmonically rich signals like speech.

**Spectrogram** of the sequence "159#" (70 ms tones, 70 ms pauses, 8 kHz).
Each digit appears as a pair of horizontal bands — one row frequency and
one column frequency.  Dashed lines mark the eight DTMF frequencies:

\htmlonly
<iframe src="dtmf_spectrogram.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

---

## ITU-T Q.24 timing constraints

[ITU-T Recommendation Q.24](https://en.wikipedia.org/wiki/Q.24) specifies
minimum timing for reliable DTMF signalling:

| Parameter                        | Minimum |
|:---------------------------------|:-------:|
| Tone duration for valid digit    | 40 ms   |
| Inter-digit pause                | 40 ms   |

In practice, telephone systems use 70--120 ms tones and pauses.
The miniDSP detector enforces the 40 ms minimums via a frame-counting
state machine; the generator asserts that requested durations meet
the minimums.

---

## Signal model

A single DTMF digit is the sum of two sinusoids at equal amplitude:

\f[
x[n] = A\,\sin\!\bigl(2\pi\, f_{\text{row}}\, n / f_s\bigr)
     + A\,\sin\!\bigl(2\pi\, f_{\text{col}}\, n / f_s\bigr),
\qquad n = 0, 1, \ldots, N_{\text{tone}}-1
\f]

where \f$A = 0.5\f$ so the peak combined amplitude is 1.0,
\f$f_s\f$ is the sampling rate, and \f$N_{\text{tone}}\f$ is the
number of samples per tone.

**Reading the formula in C:**

```c
// A -> 0.5, f_row/f_col -> row_freq/col_freq, fs -> sample_rate
// n -> i, x[n] -> output[offset + i]
for (unsigned i = 0; i < tone_samples; i++) {
    double t = (double)i / sample_rate;
    output[offset + i] = 0.5 * sin(2 * M_PI * row_freq * t)
                        + 0.5 * sin(2 * M_PI * col_freq * t);
}
```

The library implementation uses MD_sine_wave() to generate each
component separately, then sums them.

---

## Generation

**API:**

```c
unsigned len = MD_dtmf_signal_length(num_digits, sample_rate,
                                     tone_ms, pause_ms);
double *sig = malloc(len * sizeof(double));
MD_dtmf_generate(sig, "5551234", sample_rate, tone_ms, pause_ms);
```

The total signal length in samples is:

\f[
N = D \cdot \left\lfloor \frac{t_{\text{tone}} \cdot f_s}{1000} \right\rfloor
  + (D - 1) \cdot \left\lfloor \frac{t_{\text{pause}} \cdot f_s}{1000} \right\rfloor
\f]

where \f$D\f$ is the number of digits.

**Reading the formula in C:**

```c
// D -> num_digits, t_tone -> tone_ms, t_pause -> pause_ms, fs -> sample_rate
// floor(t_tone * fs / 1000) -> tone_samples
unsigned tone_samples  = (unsigned)(tone_ms * sample_rate / 1000.0);
unsigned pause_samples = (unsigned)(pause_ms * sample_rate / 1000.0);
unsigned N = num_digits * tone_samples
           + (num_digits - 1) * pause_samples;
```

**Quick example** -- generate a DTMF sequence and save as WAV:

\snippet dtmf_detector.c generate-wav

---

## Detection algorithm

Detection slides a Hanning-windowed FFT frame across the audio signal:

1. **FFT size** is the largest power of two whose window fits within
   35 ms (e.g. \f$N = 256\f$ at 8 kHz, giving \f$\Delta f = 31.25\f$ Hz).
   Keeping the window shorter than the 40 ms Q.24 minimum pause
   ensures the state machine can resolve inter-digit gaps.
2. **Hop** is \f$N/4\f$ (75 % overlap).
3. **Per frame:** apply Hanning window, compute MD_magnitude_spectrum(),
   normalise to single-sided amplitude, then check the magnitude at each
   of the eight DTMF frequency bins.
4. A digit is detected when both the strongest row and strongest column
   exceed a threshold (8\f$\times\f$ the mean spectral magnitude, roughly
   18 dB above the noise floor).
5. A **state machine** enforces ITU-T Q.24 timing:

| State       | Transition condition                    | Action                       |
|:------------|:----------------------------------------|:-----------------------------|
| **IDLE**    | Digit detected                          | Enter PENDING, start counter |
| **PENDING** | Same digit for \f$\geq\f$ 40 ms        | Enter ACTIVE (confirmed)     |
| **PENDING** | Different digit or silence              | Return to IDLE               |
| **ACTIVE**  | Same digit continues                    | Update end time              |
| **ACTIVE**  | Silence / different for \f$\geq\f$ 40 ms | Emit tone, return to IDLE   |

**Single-sided amplitude normalisation:**

\f[
\hat{X}[k] = \begin{cases}
  |X[k]| / N                   & k = 0 \text{ or } k = N/2 \\[4pt]
  2\,|X[k]| / N                & 0 < k < N/2
\end{cases}
\f]

**Reading the normalisation in C:**

```c
// X[k] -> mag[k] (raw FFTW output), N -> FFT size
for (unsigned k = 0; k < num_bins; k++) {
    mag[k] /= (double)N;                // divide by FFT size
    if (k > 0 && k < N / 2)
        mag[k] *= 2.0;                  // fold negative frequencies
}
```

**Quick example** -- detect DTMF tones in a WAV file:

\snippet dtmf_detector.c detect-file

---

## Self-test mode

Running the example with no arguments generates a known digit sequence,
detects it, and verifies correctness:

\snippet dtmf_detector.c self-test

---

## Frequency resolution and bin mapping

For a given FFT size \f$N\f$ and sampling rate \f$f_s\f$, each bin
\f$k\f$ corresponds to frequency:

\f[
f_k = k \cdot \frac{f_s}{N}
\f]

**Reading the formula in C:**

```c
// k -> bin index, fs -> sample_rate, N -> FFT size
// f_k -> freq (the frequency that bin k represents)
double freq = (double)k * sample_rate / (double)N;
```

The nearest bin for a DTMF frequency \f$f\f$ is:

\f[
k = \mathrm{round}\!\left(\frac{f \cdot N}{f_s}\right)
\f]

**Reading the formula in C:**

```c
// f -> dtmf_freq, N -> FFT size, fs -> sample_rate
// k -> bin (nearest FFT bin for the DTMF frequency)
unsigned bin = (unsigned)(dtmf_freq * N / sample_rate + 0.5);
```

The detector checks bins \f$k-1\f$, \f$k\f$, and \f$k+1\f$ and takes
the maximum magnitude, compensating for the slight frequency mismatch
when the DTMF frequency does not fall exactly on a bin centre.

At 8 kHz with \f$N = 256\f$:

| DTMF freq | Nearest bin | Bin freq  | Error   |
|:----------|:-----------:|:---------:|:-------:|
| 697 Hz    |     22      | 687.5 Hz  | -9.5    |
| 770 Hz    |     25      | 781.3 Hz  | +11.3   |
| 852 Hz    |     27      | 843.8 Hz  | -8.2    |
| 941 Hz    |     30      | 937.5 Hz  | -3.5    |
| 1209 Hz   |     39      | 1218.8 Hz | +9.8    |
| 1336 Hz   |     43      | 1343.8 Hz | +7.8    |
| 1477 Hz   |     47      | 1468.8 Hz | -8.2    |
| 1633 Hz   |     52      | 1625.0 Hz | -8.0    |

All errors are well within the ±1.5 % tolerance specified by ITU-T.
The detector also checks the two adjacent bins (±1) to handle
residual frequency mismatch.
