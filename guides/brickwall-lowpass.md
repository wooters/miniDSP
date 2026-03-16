# Brickwall Lowpass Filter {#brickwall-lowpass}

A [brickwall lowpass filter](https://en.wikipedia.org/wiki/Low-pass_filter)
removes all frequency content above a cutoff frequency with zero transition
bandwidth.  Unlike FIR or IIR filters that roll off gradually, a brickwall
filter has a perfectly sharp boundary: everything at or below the cutoff
passes through unchanged, and everything above is eliminated completely.

miniDSP implements this by transforming the signal into the frequency domain
with an FFT, zeroing the bins above the cutoff, and inverse-FFTing back
to the time domain.

---

## Frequency-domain operation

Given a real signal \f$x[n]\f$ of length \f$N\f$ with DFT \f$X(k)\f$,
the brickwall lowpass at cutoff frequency \f$f_c\f$ produces:

\f[
  X'(k) = \begin{cases}
    X(k) & \text{if } k \leq \lfloor f_c \cdot N / f_s \rfloor \\
    0    & \text{otherwise}
  \end{cases}
\f]

where \f$f_s\f$ is the sample rate.  The filtered signal is
\f$x'[n] = \mathrm{IFFT}(X'(k)) / N\f$ (the \f$1/N\f$ factor normalizes
FFTW's unnormalized inverse transform).

**Reading the formula in C:**

```c
// k -> frequency bin index, cutoff_bin -> floor(fc * N / fs)
// freq[k] -> X(k), zeroed when k > cutoff_bin
unsigned cutoff_bin = (unsigned)(cutoff_hz * N / sample_rate);
for (unsigned k = cutoff_bin + 1; k < num_bins; k++) {
    freq[k] = 0.0;   // zero both real and imaginary parts
}

// After inverse FFT, normalize every sample by N
for (unsigned n = 0; n < N; n++) {
    signal[n] /= (double)N;
}
```

---

## API

**Quick example:**

```c
double buf[4096];
MD_sine_wave(buf, 4096, 1.0, 440.0, 48000.0);
MD_lowpass_brickwall(buf, 4096, 8000.0, 48000.0);
// buf now contains only content at or below 8 kHz
```

The full signature is:

```c
void MD_lowpass_brickwall(double *signal, unsigned len,
                          double cutoff_hz, double sample_rate);
```

The filter modifies the buffer in-place.  It creates one-off FFTW plans
(`FFTW_ESTIMATE`) internally, so there is no plan caching overhead — but
it is best suited for offline/batch use rather than real-time sample-by-sample
processing.

---

## Gibbs ringing

Because the brickwall filter has an infinitely sharp cutoff in the frequency
domain, the time-domain impulse response has infinite length and exhibits
[Gibbs ringing](https://en.wikipedia.org/wiki/Gibbs_phenomenon) — oscillatory
overshoot and undershoot near the cutoff frequency.  This is the fundamental
tradeoff: zero transition bandwidth comes at the cost of time-domain ringing.

In practice, this matters only when the cutoff is close to a frequency band
you care about.  When the cutoff is far from the band of interest (e.g.,
filtering at 8 kHz when your signal content is below 4 kHz), the ringing
is inaudible and invisible.

---

## Example output

The plot below shows a 400 Hz + 12000 Hz signal before and after applying
a brickwall lowpass at 4000 Hz.  The 12000 Hz component is completely
eliminated while the 400 Hz component passes through unchanged.

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="brickwall.html" style="width:100%;height:580px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

\snippet brickwall.c apply-brickwall
