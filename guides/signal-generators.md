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
x[n] = A \sin\!\left(2\pi\!\left(f_0\,t + \tfrac{1}{2}\,\frac{f_1 - f_0}{T}\,t^2\right)\right),
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

*Future sections — Square/Sawtooth — will appear here
as new `##` headings when those generators are implemented.*
