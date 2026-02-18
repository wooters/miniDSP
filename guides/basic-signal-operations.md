# Basic Signal Operations {#basic-signal-operations}

Five fundamental time-domain operations that complement the existing signal
measurements (energy, power, entropy) and FFT-based spectrum analysis.
Together they give you a toolkit for analysing and manipulating signals
before or instead of going to the frequency domain.

---

## RMS (Root Mean Square)

RMS is the standard measure of signal "loudness" — the square root of
the mean squared amplitude:

\f[
\mathrm{RMS} = \sqrt{\frac{1}{N}\sum_{n=0}^{N-1} x[n]^2}
\f]

**Reading the formula in C:**

```c
// x[n]^2  ->  a[i] * a[i]          square each sample
// sum     ->  accumulate in a loop  sum of squares
// 1/N     ->  / (double)N           divide by number of samples
// sqrt()  ->  sqrt()                take the square root
double sum = 0.0;
for (unsigned i = 0; i < N; i++)
    sum += a[i] * a[i];
double rms = sqrt(sum / (double)N);
```

A unit-amplitude sine wave has RMS = \f$1/\sqrt{2} \approx 0.707\f$.
A DC signal of value \f$c\f$ has RMS = \f$|c|\f$.

**API:**

```c
double MD_rms(const double *a, unsigned N);
```

**Quick example** — measure the RMS of three signals:

\snippet basic_operations.c rms-measurements

**Verification tip:** `MD_rms(a, N)` should always equal
`sqrt(MD_power(a, N))` to within floating-point precision.

---

## Zero-Crossing Rate

The zero-crossing rate (ZCR) counts how often the signal changes sign,
normalised by the number of adjacent pairs:

\f[
\mathrm{ZCR} = \frac{1}{N-1}\sum_{n=1}^{N-1}
    \mathbf{1}\!\bigl[\mathrm{sgn}(x[n]) \ne \mathrm{sgn}(x[n-1])\bigr]
\f]

**Reading the formula in C:**

```c
// Walk through the signal and count sign changes.
// Zero is treated as non-negative (standard convention).
unsigned crossings = 0;
for (unsigned i = 1; i < N; i++)
    if ((a[i] < 0.0) != (a[i - 1] < 0.0))
        crossings++;
double zcr = (double)crossings / (double)(N - 1);
```

A pure sine at frequency \f$f\f$ with sample rate \f$f_s\f$ has
ZCR \f$\approx 2f/f_s\f$.  White noise has a higher ZCR than a
low-frequency tone — this makes ZCR a simple proxy for distinguishing
noise from tonal content.

**API:**

```c
double MD_zero_crossing_rate(const double *a, unsigned N);
```

**Quick example** — compare ZCR of sine, noise, and their mix:

\snippet basic_operations.c zcr-measurements

**Verification tip:** for a sine wave with an integer number of cycles
in the buffer, the ZCR should be very close to \f$2f/f_s\f$.

---

## Autocorrelation

The autocorrelation measures the similarity between a signal and a
delayed copy of itself.  The normalised form divides by \f$R[0]\f$
(the signal energy) so the output is bounded:

\f[
R[\tau] = \frac{\displaystyle\sum_{n=0}^{N-1-\tau} x[n]\,x[n+\tau]}
               {\displaystyle\sum_{n=0}^{N-1} x[n]^2}
\f]

**Reading the formula in C:**

```c
// Denominator: sum of x[n]^2  (the signal energy, R[0])
double r0 = 0.0;
for (unsigned n = 0; n < N; n++)
    r0 += a[n] * a[n];

// Numerator for each lag tau: sum of x[n] * x[n + tau]
for (unsigned tau = 0; tau < max_lag; tau++) {
    double sum = 0.0;
    for (unsigned n = 0; n < N - tau; n++)
        sum += a[n] * a[n + tau];
    out[tau] = sum / r0;   // normalise so out[0] = 1.0
}
```

\f$R[0] = 1.0\f$ always, and \f$|R[\tau]| \le 1.0\f$.
For a periodic signal, the autocorrelation peaks at the fundamental
period — this is the basis of time-domain pitch detection.

**API:**

```c
void MD_autocorrelation(const double *a, unsigned N,
                        double *out, unsigned max_lag);
```

**Quick example** — autocorrelation of a noisy sine:

\snippet basic_operations.c autocorrelation

**Verification tip:** for a 200 Hz sine at 8000 Hz sample rate, the
autocorrelation should peak at lag 40 (= 8000 / 200).

---

## Peak Detection

Peak detection finds local maxima in a signal.  A sample \f$a[i]\f$
is a peak if it is strictly greater than both immediate neighbours and
above a given threshold.  The `min_distance` parameter suppresses
nearby secondary peaks.

**Reading the algorithm:**

```c
// Left-to-right scan:
// 1. a[i] > a[i-1] AND a[i] > a[i+1]  (strict local maximum)
// 2. a[i] >= threshold                  (above noise floor)
// 3. distance from last peak >= min_distance  (suppress echoes)
```

Endpoints (i=0 and i=N-1) are never peaks because they lack two
neighbours.  A flat signal (all values equal) produces zero peaks.

**API:**

```c
void MD_peak_detect(const double *a, unsigned N, double threshold,
                    unsigned min_distance, unsigned *peaks_out,
                    unsigned *num_peaks_out);
```

**Quick example** — find peaks in the autocorrelation to estimate pitch:

\snippet basic_operations.c peak-detection

**Verification tip:** pass a hand-crafted signal like
`{0, 1, 3, 1, 0, 2, 5, 2, 0}` with threshold 0 and min_distance 1.
You should get peaks at indices 2 and 6 (values 3 and 5).

---

## Signal Mixing

Mixing computes the element-wise weighted sum of two signals:

\f[
\mathrm{out}[n] = w_a \cdot a[n] + w_b \cdot b[n]
\f]

**Reading the formula in C:**

```c
for (unsigned i = 0; i < N; i++)
    out[i] = w_a * a[i] + w_b * b[i];
```

Equal weights of 0.5 produce the average.  A weight of 1.0/0.0 passes
one signal through unchanged.  The output buffer may alias either
input (in-place safe).

**API:**

```c
void MD_mix(const double *a, const double *b, double *out,
            unsigned N, double w_a, double w_b);
```

**Quick example** — mix 80% sine with 20% noise:

\snippet basic_operations.c mix-signals

**Verification tip:** for uncorrelated signals (e.g. a sine and white
noise), the energy of the mix is approximately the sum of energies:
\f$E(\mathrm{mix}) \approx w_a^2 E(a) + w_b^2 E(b)\f$.
