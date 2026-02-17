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

*Future sections — White noise, Impulse, Chirp, Square/Sawtooth — will appear here
as new `##` headings when those generators are implemented.*
