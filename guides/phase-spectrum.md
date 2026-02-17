# Phase Spectrum {#phase-spectrum}

This tutorial walks through `examples/phase_spectrum.c`, which generates a
three-component test signal, computes its phase spectrum with
MD_phase_spectrum(), and visualises both the magnitude and phase with an
interactive Plotly chart.

If you haven't read \ref magnitude-spectrum yet, start there -- it covers
the DFT fundamentals this tutorial builds on.

## What is the phase spectrum?

Every DFT coefficient \f$X(k)\f$ is a complex number with a magnitude
and an angle.  The **magnitude spectrum** \f$|X(k)|\f$ tells you how
much energy is at frequency \f$k\f$.  The **phase spectrum**
\f$\phi(k)\f$ tells you the *timing* of that frequency component:

\f[
\phi(k) = \arg X(k) = \mathrm{atan2}\!\bigl(\mathrm{Im}\, X(k),\,\mathrm{Re}\, X(k)\bigr)
\qquad \phi(k) \in [-\pi,\, \pi]
\f]

### Intuitive meaning

| Signal | Phase at its bin |
|--------|-----------------|
| \f$\cos(2\pi k_0 n / N)\f$ | \f$\phi(k_0) = 0\f$ |
| \f$\sin(2\pi k_0 n / N)\f$ | \f$\phi(k_0) = -\pi/2\f$ |
| \f$\cos(2\pi k_0 n / N + \varphi)\f$ | \f$\phi(k_0) = \varphi\f$ |
| Impulse delayed by \f$d\f$ samples | \f$\phi(k) = -2\pi k d / N\f$ (linear) |

A **time-delayed signal** is the most important case: the DFT shift
theorem says delaying a signal by \f$d\f$ samples adds a linear phase
ramp of slope \f$-2\pi d / N\f$ radians per bin.  This is the basis
of the GCC-PHAT delay estimator (see MD_get_delay()).

### When to trust the phase

Phase is only meaningful at bins where the magnitude is significant.
At bins dominated by noise or leakage, \f$\phi(k)\f$ is numerically
unreliable -- always examine MD_magnitude_spectrum() alongside the
phase spectrum.

## Step 1: Generate a test signal {#phase-spectrum-step1}

The test signal has three components at **exact integer-bin frequencies**
(N = 1024, sample rate = 1024 Hz, so bin \f$k\f$ corresponds to exactly
\f$k\f$ Hz).  With exact bins there is no spectral leakage, so no
windowing is needed and the phase values are bit-exact.

\snippet phase_spectrum.c generate-signal

The three components and their expected phase values:

| Component | Bin | Expected \f$\phi\f$ |
|-----------|-----|---------------------|
| `cos(50 Hz)` | 50 | 0 |
| `sin(100 Hz)` | 100 | \f$-\pi/2\f$ |
| `cos(200 Hz, \f$\pi/4\f$ offset)` | 200 | \f$\pi/4\f$ |

## Step 2: Compute the phase spectrum {#phase-spectrum-step2}

MD_phase_spectrum() reuses the same cached FFT plan as
MD_magnitude_spectrum() and MD_power_spectral_density().
Calling both in sequence costs only one plan lookup:

\snippet phase_spectrum.c compute-phase

The output `phase[k]` is in radians, range \f$[-\pi, \pi]\f$.  The
magnitude output `mag[k]` lets you distinguish signal bins from
noise bins when interpreting the phase.

## Results {#phase-spectrum-results}

\image html phase-spectrum.png "Phase spectrum (top: magnitude, bottom: phase)" width=700px

The magnitude plot (top) confirms that energy is concentrated at bins 50,
100, and 200.  The phase plot (bottom) shows:

- Bin 50: \f$\phi \approx 0\f$ (pure cosine)
- Bin 100: \f$\phi \approx -\pi/2\f$ (pure sine)
- Bin 200: \f$\phi \approx \pi/4\f$ (phase-shifted cosine)

At all other bins the phase is numerically arbitrary (magnitude is zero).

## Key takeaways {#phase-spectrum-takeaways}

- \f$\phi(k) = \arg X(k)\f$ is computed by `carg()`, which calls
  `atan2(Im, Re)` and returns values in \f$[-\pi, \pi]\f$.
- Phase is **scale-invariant**: multiplying a signal by a positive
  constant does not change its phase.
- **No windowing** is needed when all frequency components land on
  exact integer bins.  Windowing smears phase and should be omitted
  for clean phase measurements.
- Always look at the magnitude spectrum alongside the phase spectrum.
  Phase at low-magnitude bins is meaningless noise.
- A linear phase ramp (\f$\phi(k) = -2\pi k d/N\f$) indicates a
  time delay of \f$d\f$ samples -- the foundation of GCC-PHAT delay
  estimation.

## Connection to phase vocoder

The phase vocoder uses phase differences between consecutive STFT
frames to track the instantaneous frequency of each bin, enabling
high-quality time-stretching and pitch-shifting of audio.
MD_phase_spectrum() provides the per-frame phase needed as input
to such an algorithm.

## Further reading {#phase-spectrum-reading}

- [Phase (waves)](https://en.wikipedia.org/wiki/Phase_(waves)) -- Wikipedia
- [DFT shift theorem](https://en.wikipedia.org/wiki/Discrete_Fourier_transform#Shift_theorem) -- Wikipedia
- [Phase vocoder](https://en.wikipedia.org/wiki/Phase_vocoder) -- Wikipedia
- Julius O. Smith, [Spectral Audio Signal Processing](https://ccrma.stanford.edu/~jos/sasp/) -- free online textbook
- Previous tutorial: \ref power-spectral-density
- Next: \ref stft-spectrogram

## API reference {#phase-spectrum-api}

- MD_phase_spectrum() -- compute the one-sided phase spectrum
- MD_magnitude_spectrum() -- compute the magnitude spectrum
- MD_power_spectral_density() -- compute the PSD (periodogram)
- MD_shutdown() -- free cached FFTW plans
