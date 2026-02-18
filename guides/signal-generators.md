# Signal Generators {#signal-generators}

Signal generators produce test signals from scratch — no audio file or microphone required.
They are the "hello world" of DSP: a pure sine wave, an impulse, or a burst of noise gives you
a known input that you can trace through any processing chain and verify at every step.

miniDSP ships generators as simple, stateless functions.  They take an output buffer
and a handful of parameters; no setup or teardown is needed.

---

## Sine wave

A sine wave at frequency \f$f\f$ Hz, amplitude \f$A\f$, and sample rate \f$f_s\f$:

\f[
x[n] = A \sin\!\left(\frac{2\pi f n}{f_s}\right), \quad n = 0, 1, \ldots, N-1
\f]

**API:**

```c
void MD_sine_wave(double *output, unsigned N, double amplitude,
                  double freq, double sample_rate);
```

**Quick example** — generate one second of A4 (440 Hz):

\snippet sine_wave.c generate-signal

**Verifying via spectrum**

The clearest way to confirm the generator is correct is to feed its output to
`MD_magnitude_spectrum()` and check that the peak lands on the expected bin:

```c
unsigned N = 1024;
double   fs = 16000.0, freq = 1000.0;   /* bin 64 */
double   sig[N], mag[N/2 + 1];

MD_sine_wave(sig, N, 1.0, freq, fs);
MD_magnitude_spectrum(sig, N, mag);
/* mag[64] should be the largest value */
```

See `examples/sine_wave.c` for a full runnable program that generates the spectrum
and writes an interactive HTML plot.

---

## Impulse

A discrete impulse (Kronecker delta) at position \f$n_0\f$ with amplitude \f$A\f$:

\f[
x[n] = \begin{cases} A & \text{if } n = n_0 \\ 0 & \text{otherwise} \end{cases}
\f]

The unit impulse (\f$A = 1\f$, \f$n_0 = 0\f$) is the identity element of
convolution and has a perfectly flat magnitude spectrum — every frequency
bin equals 1.0.

**API:**

```c
void MD_impulse(double *output, unsigned N, double amplitude, unsigned position);
```

**Quick example** — generate a unit impulse at sample 0:

\snippet impulse.c generate-impulse

**Verifying via spectrum**

Feed the impulse to `MD_magnitude_spectrum()` and confirm every bin has the
same magnitude:

```c
unsigned N = 1024;
double sig[N], mag[N/2 + 1];

MD_impulse(sig, N, 1.0, 0);
MD_magnitude_spectrum(sig, N, mag);
/* Every element of mag[] should be 1.0 */
```

See `examples/impulse.c` for a full runnable program that generates both
time-domain and frequency-domain plots.

---

## Chirp (swept sine)

A **chirp** sweeps frequency over time — either linearly or logarithmically.
Chirps are the standard test signal for spectrograms and for measuring filter
magnitude response.

### Linear chirp

A linear chirp sweeps instantaneous frequency from \f$f_0\f$ to \f$f_1\f$ at a
constant rate over duration \f$T = (N-1)/f_s\f$:

\f[
x[n] = A \sin\!\left(2\pi\!\left(f_0\,t + \frac{1}{2}\,\frac{f_1 - f_0}{T}\,t^2\right)\right),
\quad t = n / f_s
\f]

**API:**

```c
void MD_chirp_linear(double *output, unsigned N, double amplitude,
                     double f_start, double f_end, double sample_rate);
```

### Logarithmic chirp

A logarithmic chirp sweeps frequency exponentially, spending equal time per
octave.  This is ideal for audio systems whose response is best viewed on a
log-frequency axis.

\f[
x[n] = A \sin\!\left(\frac{2\pi f_0 T}{\ln r}\!\left(r^{t/T} - 1\right)\right),
\quad r = f_1 / f_0,\quad t = n / f_s
\f]

**API:**

```c
void MD_chirp_log(double *output, unsigned N, double amplitude,
                  double f_start, double f_end, double sample_rate);
```

Requires \f$f_0 > 0\f$, \f$f_1 > 0\f$, and \f$f_0 \ne f_1\f$.

**Quick example** — generate both chirp types and compare:

\snippet chirp_wave.c generate-chirps

See `examples/chirp_wave.c` for a full runnable program that generates the
magnitude spectra of both chirp types and writes an interactive HTML plot.

---

## Square wave

A square wave at frequency \f$f\f$ Hz alternates between \f$+A\f$ and \f$-A\f$:

\f[
x[n] = \begin{cases}
  +A  & 0 < \phi < \pi \\
  -A  & \pi < \phi < 2\pi \\
   0  & \phi = 0 \text{ or } \phi = \pi
\end{cases}
\f]

where \f$\phi = 2\pi f n / f_s \pmod{2\pi}\f$.

Its Fourier series contains **only odd harmonics** (1f, 3f, 5f, …) with amplitudes
decaying as \f$1/k\f$ — a textbook demonstration of the Gibbs phenomenon.

**API:**

```c
void MD_square_wave(double *output, unsigned N, double amplitude,
                    double freq, double sample_rate);
```

**Quick example:**

\snippet square_sawtooth.c generate-square

See `examples/square_sawtooth.c` for a full program that compares the square and
sawtooth spectra side by side.

---

## Sawtooth wave

A sawtooth wave ramps linearly from \f$-A\f$ to \f$+A\f$ over each period:

\f[
x[n] = A \left(\frac{\phi}{\pi} - 1\right)
\f]

where \f$\phi = 2\pi f n / f_s \pmod{2\pi}\f$.

Unlike the square wave, the sawtooth contains **all integer harmonics**
(1f, 2f, 3f, …), each decaying as \f$1/k\f$.  Comparing the two spectra
shows how waveform shape determines harmonic content.

**API:**

```c
void MD_sawtooth_wave(double *output, unsigned N, double amplitude,
                      double freq, double sample_rate);
```

**Quick example:**

\snippet square_sawtooth.c generate-sawtooth

See `examples/square_sawtooth.c` for the full runnable program.

---

## White noise

Gaussian white noise has **equal power at all frequencies** — its power spectral
density (PSD) is approximately flat across the entire band.  It is the standard
broadband test signal for filter characterisation, impulse response measurement,
and SNR experiments.

Each sample is drawn independently from a normal distribution with mean 0 and
standard deviation \f$\sigma\f$:

\f[
x[n] \sim \mathrm{N}(0,\, \sigma^2), \quad n = 0, 1, \ldots, N-1
\f]

Samples are generated with the [Box-Muller transform](https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform) seeded by `seed`, so the
same seed always produces the same sequence — useful for reproducible tests.

**API:**

```c
void MD_white_noise(double *output, unsigned N, double amplitude,
                    unsigned seed);
```

`amplitude` is the standard deviation \f$\sigma\f$ of the distribution.

**Quick example** — generate 4096 samples of unit-variance noise:

\snippet white_noise.c generate-signal

**Verifying via PSD**

Feed the noise to `MD_power_spectral_density()` and confirm the spectrum is
approximately flat — no bin should dominate:

```c
unsigned N = 4096;
double *sig = malloc(N * sizeof(double));
double *psd = malloc((N/2 + 1) * sizeof(double));

MD_white_noise(sig, N, 1.0, 42);
MD_power_spectral_density(sig, N, psd);
/* psd[] should fluctuate around a constant level */
```

See `examples/white_noise.c` for a full runnable program that computes and
plots the power spectral density.
