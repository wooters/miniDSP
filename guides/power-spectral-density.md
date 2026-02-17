# Computing the Power Spectral Density {#power-spectral-density}

This tutorial walks through the `examples/power_spectral_density.c` program,
which generates a multi-tone signal, windows it, computes the PSD with
MD_power_spectral_density(), and visualises the result.

If you haven't already, read the \ref magnitude-spectrum tutorial first --
it covers the DFT fundamentals that this tutorial builds on.

## What is the power spectral density?

The **magnitude spectrum** (covered in the \ref magnitude-spectrum tutorial)
shows *amplitude* at each frequency. The **power spectral density (PSD)**
shows how the signal's *power* is distributed across frequencies.

Given the DFT coefficients \f$X(k)\f$ of an \f$N\f$-sample signal, the
periodogram estimate of the PSD is:

\f[
\text{PSD}(k) = \frac{|X(k)|^2}{N}
\f]

The key difference from the magnitude spectrum is the **squaring** --
we use \f$|X(k)|^2\f$ (power) instead of \f$|X(k)|\f$ (amplitude).

## Magnitude spectrum vs. PSD

| Property | Magnitude spectrum | PSD |
|----------|-------------------|-----|
| Formula  | \f$\|X(k)\|\f$ | \f$\|X(k)\|^2 / N\f$ |
| Units    | Amplitude | Power |
| dB conversion | \f$20 \log_{10}\f$ | \f$10 \log_{10}\f$ |
| Use case | Identify frequency components | Noise analysis, SNR estimation |

Why the different dB formulas? Power is proportional to the *square*
of amplitude. Since \f$10 \log_{10}(A^2) = 20 \log_{10}(A)\f$, using
\f$10 \log_{10}\f$ for power and \f$20 \log_{10}\f$ for amplitude
produces the same dB values for equivalent signals.

## Parseval's theorem

Parseval's theorem states that the total energy in the time domain equals
the total energy in the frequency domain:

\f[
\sum_{n=0}^{N-1} |x[n]|^2 = \frac{1}{N} \sum_{k=0}^{N-1} |X(k)|^2
  = \sum_{k=0}^{N-1} \text{PSD}(k)
\f]

This is a useful sanity check: the sum of all PSD bins should equal
\f$\sum x[n]^2\f$ (the signal's total energy). The miniDSP test suite
verifies this property.

## Step 1: Generate a test signal

The test signal is identical to the magnitude spectrum example: three
sinusoids at 440 Hz, 1000 Hz, and 2500 Hz, plus a DC offset.

\snippet power_spectral_density.c generate-signal

## Step 2: Apply a window function

Windowing is just as important for PSD estimation as for the magnitude
spectrum. Without it, spectral leakage smears power across bins,
making it hard to distinguish signal from noise.

\snippet power_spectral_density.c apply-window

## Step 3: Compute the PSD

A single call to MD_power_spectral_density() computes
\f$\text{PSD}(k) = |X(k)|^2 / N\f$ for bins \f$k = 0, 1, \ldots, N/2\f$.
The \f$/N\f$ normalisation is handled internally.

\snippet power_spectral_density.c compute-psd

## Step 4: Convert to a one-sided PSD

Just like the magnitude spectrum, we convert to a one-sided
representation by doubling the interior bins. DC and Nyquist bins
are not doubled since they have no mirror image.

Note that we do **not** divide by \f$N\f$ again -- MD_power_spectral_density()
already includes the \f$/N\f$ normalisation.

\snippet power_spectral_density.c one-sided-psd

## Results

### Linear scale

The linear PSD plot shows power concentrated at the three input
frequencies. Because PSD squares the amplitude, the ratio between
peaks is exaggerated compared to the magnitude spectrum: the 440 Hz
tone (amplitude 1.0, power 1.0) towers over the 2500 Hz tone
(amplitude 0.3, power 0.09).

\image html power-spectral-density-linear.png "Power spectral density -- linear scale" width=700px

### Logarithmic (dB) scale

The dB plot uses \f$10 \log_{10}(\text{PSD})\f$ -- note the factor of
10 (not 20) because PSD is already a power quantity. This reveals the
noise floor and window side lobes clearly.

\image html power-spectral-density-db.png "Power spectral density -- dB scale" width=700px

## Key takeaways

- PSD measures how **power** distributes across frequencies; the
  magnitude spectrum measures **amplitude**.
- Use \f$10 \log_{10}\f$ for PSD (power) and \f$20 \log_{10}\f$ for
  magnitude spectrum (amplitude) -- both give the same dB values for
  equivalent signals.
- **Parseval's theorem** provides a sanity check: the sum of all PSD
  bins equals the signal's total energy \f$\sum x[n]^2\f$.
- PSD is widely used in noise analysis, SNR estimation, and anywhere
  you need to quantify the power distribution of a signal.

## Further reading

- [Spectral density](https://en.wikipedia.org/wiki/Spectral_density) -- Wikipedia
- [Parseval's theorem](https://en.wikipedia.org/wiki/Parseval%27s_theorem) -- Wikipedia
- Julius O. Smith, [Spectral Audio Signal Processing](https://ccrma.stanford.edu/~jos/sasp/) -- free online textbook
- Previous tutorial: \ref magnitude-spectrum

## API reference

- MD_power_spectral_density() -- compute the PSD (periodogram)
- MD_magnitude_spectrum() -- compute the magnitude spectrum
- MD_Gen_Hann_Win() -- generate a Hanning window
- MD_shutdown() -- free cached FFTW plans
