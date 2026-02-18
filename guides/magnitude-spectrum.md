# Computing the Magnitude Spectrum {#magnitude-spectrum}

This tutorial walks through the `examples/magnitude_spectrum.c` program,
which generates a multi-tone signal, windows it, computes the magnitude
spectrum with MD_magnitude_spectrum(), and visualises the result.

## What is a magnitude spectrum?

Every real-valued signal can be decomposed into a sum of sinusoids at
different frequencies. The **magnitude spectrum** tells you the amplitude
of each sinusoidal component.

The tool for this decomposition is the **Discrete Fourier Transform (DFT)**.
Given \f$N\f$ samples \f$x[n]\f$, the DFT produces \f$N\f$ complex
coefficients:

\f[
X(k) = \sum_{n=0}^{N-1} x[n] \, e^{-j\,2\pi\,k\,n/N}, \qquad k = 0, 1, \ldots, N-1
\f]

The **magnitude spectrum** is simply the absolute value of these
coefficients: \f$|X(k)|\f$. Each bin \f$k\f$ corresponds to frequency
\f$f_k = k \cdot f_s / N\f$, where \f$f_s\f$ is the sample rate.

For a real signal, the spectrum is symmetric around \f$N/2\f$, so only
the first \f$N/2 + 1\f$ bins (the "one-sided" spectrum) carry unique
information.

## Step 1: Generate a test signal

We create a signal with three known sinusoidal components (440 Hz,
1000 Hz, and 2500 Hz) at different amplitudes, plus a small DC offset.
This lets us verify that the spectrum shows peaks at exactly those
frequencies.

\snippet magnitude_spectrum.c generate-signal

## Step 2: Apply a window function

Before computing the FFT, we multiply the signal by a **[Hanning window](https://en.wikipedia.org/wiki/Hann_function)**.
Why?

The DFT assumes the input is one period of a periodic signal. In reality,
our signal chunk rarely starts and ends at a perfect zero crossing. This
mismatch creates an artificial discontinuity at the boundaries, which
causes energy to "leak" from the true frequency into neighbouring bins --
an effect called **spectral leakage**.

A window function tapers the signal smoothly to zero at both edges,
eliminating the discontinuity:

\f[
w[n] = 0.5 \left(1 - \cos\!\left(\frac{2\pi\,n}{N-1}\right)\right)
\f]

\snippet magnitude_spectrum.c apply-window

The trade-off: windowing widens the main lobe of each peak slightly
(reducing frequency resolution) but dramatically suppresses the side
lobes (improving dynamic range).

## Step 3: Compute the magnitude spectrum

With the windowed signal ready, a single call to MD_magnitude_spectrum()
computes \f$|X(k)|\f$ for bins \f$k = 0, 1, \ldots, N/2\f$:

\snippet magnitude_spectrum.c compute-magnitude

Internally, miniDSP uses FFTW to compute the real-to-complex FFT, then
takes the absolute value of each complex coefficient.

## Step 4: Convert to a one-sided amplitude spectrum

The raw magnitudes \f$|X(k)|\f$ are unnormalised (they scale with \f$N\f$).
To recover the actual signal amplitudes, we:

1. **Divide every bin by \f$N\f$** to undo the DFT scaling.
2. **Double the interior bins** (\f$k = 1 \ldots N/2-1\f$) because the
   energy from the discarded negative-frequency bins folds onto the
   positive side.
3. **Leave DC and Nyquist alone** -- they have no mirror image.

\snippet magnitude_spectrum.c one-sided-amplitude

After this normalisation, a pure sine wave at amplitude \f$A\f$ produces
a peak of height \f$A\f$ in the spectrum.

## Results

### Linear scale

The linear-scale plot shows clear spikes at the three input frequencies.
The 440 Hz tone (amplitude 1.0) is tallest, followed by 1000 Hz (0.6)
and 2500 Hz (0.3). The tiny spike at 0 Hz is the DC offset (0.1).

\image html magnitude-spectrum-linear.png "Magnitude spectrum -- linear scale" width=700px

### Logarithmic (dB) scale

The dB plot reveals low-level details invisible on the linear scale.
The Hanning window's side lobes appear as the gradually decaying skirt
around each peak. Without windowing, these lobes would be much higher
and wider, obscuring nearby weak signals.

\image html magnitude-spectrum-db.png "Magnitude spectrum -- dB scale" width=700px

## Key takeaways

- The DFT converts a time-domain signal into its frequency components.
- Always **window** the signal before the FFT to control spectral leakage.
- **Normalise** the output by dividing by \f$N\f$ and doubling interior
  bins to get a one-sided amplitude spectrum.
- Use a **dB scale** (\f$20 \log_{10}\f$) to see low-level details.

## Further reading

- [Discrete Fourier transform](https://en.wikipedia.org/wiki/Discrete_Fourier_transform) -- Wikipedia
- [Spectral leakage](https://en.wikipedia.org/wiki/Spectral_leakage) -- Wikipedia
- Julius O. Smith, [Mathematics of the DFT](https://ccrma.stanford.edu/~jos/mdft/) -- free online textbook
- Next tutorial: \ref power-spectral-density

## API reference

- MD_magnitude_spectrum() -- compute the magnitude spectrum
- MD_Gen_Hann_Win() -- generate a Hanning window
- MD_shutdown() -- free cached FFTW plans
