# Sample Rate Conversion {#resampling}

[Sample rate conversion](https://en.wikipedia.org/wiki/Sample-rate_conversion) (resampling)
changes a signal from one sample rate to another. This is essential when
combining audio from different sources, preparing data for systems that
expect a specific rate, or reducing storage by downsampling.

miniDSP provides a high-quality offline [polyphase](https://en.wikipedia.org/wiki/Polyphase_filter)
sinc resampler that handles arbitrary rate ratios (e.g., 44100 Hz to 48000 Hz)
with >100 dB stopband attenuation using default parameters.

---

## Mathematical building blocks

### Bessel I₀

The zeroth-order modified Bessel function of the first kind appears in the
[Kaiser window](https://en.wikipedia.org/wiki/Kaiser_window) formula.
It is computed via the convergent power series:

\f[
I_0(x) = \sum_{k=0}^{\infty} \left[\frac{(x/2)^k}{k!}\right]^2
\f]

**Reading the formula in C:**

```c
// x -> input value, sum -> I₀(x), term -> (x/2)^k / k!
double sum = 1.0;
double term = 1.0;
double half_x = x / 2.0;

for (unsigned k = 1; k < 300; k++) {
    term *= (half_x / (double)k);
    double term_sq = term * term;
    sum += term_sq;
    if (term_sq < 1e-15 * sum) break;  // converged
}
// sum now holds I₀(x)
```

**API:**

```c
double MD_bessel_i0(double x);
```

### Normalized sinc

The [sinc function](https://en.wikipedia.org/wiki/Sinc_function) is the ideal
lowpass interpolation kernel. The normalized form is:

\f[
\mathrm{sinc}(x) = \begin{cases}
  1 & \text{if } |x| < 10^{-12} \\
  \dfrac{\sin(\pi x)}{\pi x} & \text{otherwise}
\end{cases}
\f]

**Reading the formula in C:**

```c
// x -> input value, result -> sinc(x)
double result;
if (fabs(x) < 1e-12) {
    result = 1.0;
} else {
    double px = M_PI * x;
    result = sin(px) / px;
}
```

**API:**

```c
double MD_sinc(double x);
```

---

## The polyphase sinc resampler

### How it works

The resampler converts a signal from sample rate \f$f_{\mathrm{in}}\f$ to
\f$f_{\mathrm{out}}\f$ by treating each output sample as a fractional-position
lookup into the input signal, filtered through a windowed sinc kernel.

For each output sample \f$n\f$:

1. Compute the fractional input position:
   \f$p = n \cdot f_{\mathrm{in}} / f_{\mathrm{out}}\f$
2. Split into integer index \f$\lfloor p \rfloor\f$ and fractional offset
3. Select two adjacent filter sub-phases from a precomputed table
4. Linearly interpolate the sub-phase coefficients
5. Dot product with surrounding input samples

The filter table contains 512 sub-phases, each with
\f$2 \times \mathrm{num\_zero\_crossings}\f$ taps of a Kaiser-windowed sinc.
Anti-aliasing is handled automatically: for downsampling, the sinc cutoff is
scaled to \f$\min(f_{\mathrm{in}}, f_{\mathrm{out}})/2\f$.

### Output buffer sizing

\f[
N_{\mathrm{out}} = \left\lceil N_{\mathrm{in}} \cdot \frac{f_{\mathrm{out}}}{f_{\mathrm{in}}} \right\rceil
\f]

**Reading the formula in C:**

```c
// input_len -> N_in, in_rate -> f_in, out_rate -> f_out
unsigned out_len = (unsigned)ceil((double)input_len * out_rate / in_rate);
```

**API:**

```c
unsigned MD_resample_output_len(unsigned input_len,
                                double in_rate, double out_rate);
```

### Resampling a signal

**API:**

```c
unsigned MD_resample(const double *input, unsigned input_len,
                     double *output, unsigned max_output_len,
                     double in_rate, double out_rate,
                     unsigned num_zero_crossings, double kaiser_beta);
```

**Parameters:**

- `num_zero_crossings` controls filter quality. Range: 8 (fast) to 64 (high quality). Default recommendation: 32.
- `kaiser_beta` controls stopband attenuation. Recommendation: 10.0 for ~100 dB.
- Returns the number of samples written to `output`.

**Quick example:**

\snippet resampler.c resample-basic

---

## Common rate pairs

| Conversion | Ratio | Use case |
|---|---|---|
| 44100 to 48000 | 160/147 | CD audio to professional video/DAW |
| 48000 to 44100 | 147/160 | Professional to CD quality |
| 48000 to 16000 | 1/3 | Wideband to narrowband speech |
| 44100 to 22050 | 1/2 | Half-rate for reduced storage |
| 16000 to 8000 | 1/2 | Wideband to telephone bandwidth |

---

## Further reading

- [Sample-rate conversion](https://en.wikipedia.org/wiki/Sample-rate_conversion) -- Wikipedia
- [Sinc filter](https://en.wikipedia.org/wiki/Sinc_filter) -- Wikipedia
- [Polyphase filter](https://en.wikipedia.org/wiki/Polyphase_filter) -- Wikipedia
- [Kaiser window](https://en.wikipedia.org/wiki/Kaiser_window) -- Wikipedia

## API reference

- MD_bessel_i0() -- zeroth-order modified Bessel function I₀
- MD_sinc() -- normalized sinc function
- MD_resample_output_len() -- compute required output buffer size
- MD_resample() -- polyphase sinc resampler
