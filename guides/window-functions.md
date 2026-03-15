# Window Functions {#window-functions}

[Window functions](https://en.wikipedia.org/wiki/Window_function) taper a finite signal block before an [FFT](https://en.wikipedia.org/wiki/Fast_Fourier_transform) so the block
edges do not create a large discontinuity. That discontinuity causes
**[spectral leakage](https://en.wikipedia.org/wiki/Spectral_leakage)**: energy spreads into neighboring bins.

miniDSP provides five windows so you can compare the trade-off
between main-lobe width (frequency resolution) and [sidelobe](https://en.wikipedia.org/wiki/Sidelobes) level
(leakage suppression).

---

## Hanning window

The Hanning (Hann) window is a common default:

\f[
w[n] = 0.5 \left(1 - \cos\!\left(\frac{2\pi n}{N-1}\right)\right),
\quad n = 0, 1, \ldots, N-1
\f]

It tapers smoothly to zero at both ends and gives good all-around
performance for FFT analysis.

**Reading the formula in C:**

```c
// n -> N (window length), i -> n (sample index), out[i] -> w[n]
if (n == 1) {
    out[0] = 1.0;
} else {
    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / n_minus_1));
    }
}
```

**API:**

```c
void MD_Gen_Hann_Win(double *out, unsigned n);
```

**Visuals** — window taps and magnitude response:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="hann_window_time.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="hann_window_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

The smooth taper to zero at both ends reduces leakage compared with a
rectangular window.

**Quick example:**

\snippet window_functions.c generate-hann

---

## Hamming window

The Hamming window keeps a similar shape to Hanning, but with non-zero
endpoints and a lower first sidelobe:

\f[
w[n] = 0.54 - 0.46 \cos\!\left(\frac{2\pi n}{N-1}\right)
\f]

**Reading the formula in C:**

```c
// n -> N (window length), i -> n (sample index), out[i] -> w[n]
if (n == 1) {
    out[0] = 1.0;
} else {
    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)i / n_minus_1);
    }
}
```

**API:**

```c
void MD_Gen_Hamming_Win(double *out, unsigned n);
```

**Visuals** — window taps and magnitude response:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="hamming_window_time.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="hamming_window_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

Compared with Hanning, the non-zero endpoints and coefficients shift
the sidelobe pattern while keeping a similar main-lobe width.

**Quick example:**

\snippet window_functions.c generate-hamming

---

## Blackman window

The Blackman window strongly suppresses sidelobes by adding another
cosine term:

\f[
w[n] = 0.42
     - 0.5 \cos\!\left(\frac{2\pi n}{N-1}\right)
     + 0.08 \cos\!\left(\frac{4\pi n}{N-1}\right)
\f]

Compared with Hanning/Hamming, it has much lower sidelobes but a wider
main lobe.

**Reading the formula in C:**

```c
// n -> N (window length), i -> n (sample index),
// p -> 2*pi*n/(N-1), out[i] -> w[n]
if (n == 1) {
    out[0] = 1.0;
} else {
    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        double p = 2.0 * M_PI * (double)i / n_minus_1;
        out[i] = 0.42 - 0.5 * cos(p) + 0.08 * cos(2.0 * p);
    }
}
```

**API:**

```c
void MD_Gen_Blackman_Win(double *out, unsigned n);
```

**Visuals** — window taps and magnitude response:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="blackman_window_time.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="blackman_window_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

You should see much lower sidelobes than Hanning/Hamming, with a wider
main lobe in the response plot.

**Quick example:**

\snippet window_functions.c generate-blackman

---

## Rectangular window

The rectangular window is the no-taper baseline:

\f[
w[n] = 1
\f]

It preserves the narrowest main lobe but has the highest sidelobes.

**Reading the formula in C:**

```c
// i -> n (sample index), out[i] -> w[n]
for (unsigned i = 0; i < n; i++) {
    out[i] = 1.0;
}
```

**API:**

```c
void MD_Gen_Rect_Win(double *out, unsigned n);
```

**Visuals** — window taps and magnitude response:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="rect_window_time.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="rect_window_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

As the no-taper baseline, rectangular gives the narrowest main lobe and
the highest sidelobes.

**Quick example:**

\snippet window_functions.c generate-rect

---

## Kaiser window

The [Kaiser window](https://en.wikipedia.org/wiki/Kaiser_window) uses the
zeroth-order modified Bessel function \f$I_0\f$ to provide continuous control
over the sidelobe/mainlobe tradeoff via a single parameter \f$\beta\f$:

\f[
w[n] = \frac{I_0\!\left(\beta\,\sqrt{1 - \left(\frac{2n}{N-1}-1\right)^2}\right)}
            {I_0(\beta)},
\quad n = 0, 1, \ldots, N-1
\f]

Higher \f$\beta\f$ values produce lower sidelobes (better leakage suppression)
at the cost of a wider main lobe:
- \f$\beta \approx 5\f$ : ~45 dB stopband attenuation
- \f$\beta \approx 10\f$ : ~100 dB stopband attenuation
- \f$\beta \approx 14\f$ : ~120 dB stopband attenuation

**Reading the formula in C:**

```c
// n -> N (window length), i -> n (sample index), out[i] -> w[n]
// beta -> shape parameter, inv_i0_beta -> 1/I₀(β)
if (n == 1) {
    out[0] = 1.0;
} else {
    double inv_i0_beta = 1.0 / MD_bessel_i0(beta);
    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        double t = 2.0 * (double)i / n_minus_1 - 1.0;
        double arg = 1.0 - t * t;
        if (arg < 0.0) arg = 0.0;
        out[i] = MD_bessel_i0(beta * sqrt(arg)) * inv_i0_beta;
    }
}
```

**API:**

```c
void MD_Gen_Kaiser_Win(double *out, unsigned n, double beta);
```

**Quick example:**

\snippet window_functions.c generate-kaiser

---

## Quick comparison

| Window | Edge values | Sidelobes | Main lobe | Tunable? |
|--------|-------------|-----------|-----------|----------|
| Rectangular | 1.0 | Highest | Narrowest | No |
| Hanning | 0.0 | Low | Medium | No |
| Hamming | 0.08 | Lower first sidelobe | Medium | No |
| Blackman | 0.0 | Very low | Widest | No |
| Kaiser | > 0 | Configurable via \f$\beta\f$ | Configurable | Yes |

If you are unsure where to start, Hanning is a good default. Use
Blackman when leakage suppression matters more than peak sharpness.
Use Kaiser when you need precise control over the sidelobe level
(e.g., for FIR filter design or resampling).
All response plots above use the same tap length and zero-padded FFT
size, so sidelobe and main-lobe differences are directly comparable.

## Further reading

- [Window function](https://en.wikipedia.org/wiki/Window_function) -- Wikipedia
- [Hann function](https://en.wikipedia.org/wiki/Hann_function) -- Wikipedia
- [Hamming window](https://en.wikipedia.org/wiki/Window_function#Hann_and_Hamming_windows) -- Wikipedia
- [Blackman window](https://en.wikipedia.org/wiki/Window_function#Blackman_window) -- Wikipedia
- [Kaiser window](https://en.wikipedia.org/wiki/Kaiser_window) -- Wikipedia

## API reference

- MD_Gen_Hann_Win() -- generate a Hanning window
- MD_Gen_Hamming_Win() -- generate a Hamming window
- MD_Gen_Blackman_Win() -- generate a Blackman window
- MD_Gen_Rect_Win() -- generate a rectangular window
- MD_Gen_Kaiser_Win() -- generate a Kaiser window with configurable \f$\beta\f$
