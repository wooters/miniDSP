# FIR Filters and Convolution {#fir-convolution}

[Finite impulse response (FIR)](https://en.wikipedia.org/wiki/Finite_impulse_response) filtering and [convolution](https://en.wikipedia.org/wiki/Convolution) are core DSP tools:
you shape a signal by summing delayed, weighted copies of it.

miniDSP provides four related APIs so students can start with direct
time-domain sums and then compare against the FFT overlap-add fast method.

---

## Time-domain convolution

For input \f$x[n]\f$ (length \f$N\f$) and kernel \f$h[k]\f$ (length \f$M\f$),
full linear convolution is:

\f[
y[n] = \sum_{k=0}^{M-1} h[k]\,x[n-k], \quad n = 0,1,\ldots,N+M-2
\f]

Out-of-range samples of \f$x[\cdot]\f$ are treated as zero.

**Reading the formula in C:**

```c
// n -> output index, k -> kernel index, x[n-k] -> signal[si], h[k] -> kernel[k]
for (unsigned n = 0; n < out_len; n++) {
    double acc = 0.0;
    for (unsigned k = 0; k < kernel_len; k++) {
        if (k > n) break;            // n-k would be negative
        unsigned si = n - k;
        if (si >= signal_len) continue;
        acc += signal[si] * kernel[k];
    }
    out[n] = acc;
}
```

**API:**

```c
unsigned MD_convolution_num_samples(unsigned signal_len, unsigned kernel_len);

void MD_convolution_time(const double *signal, unsigned signal_len,
                         const double *kernel, unsigned kernel_len,
                         double *out);
```

**Visuals** — response and magnitude spectrum:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="conv_time_response.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="conv_time_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

**Quick example:**

\snippet fir_convolution.c convolution-time

---

## Moving-average FIR filter

The moving-average filter is the simplest FIR [low-pass](https://en.wikipedia.org/wiki/Low-pass_filter):

\f[
y[n] = \frac{1}{M}\sum_{k=0}^{M-1} x[n-k]
\f]

This implementation is causal and uses zero-padded startup (fixed divisor
\f$M\f$ at every sample).

**Reading the formula in C:**

```c
// M -> window_len, x[n] -> signal[n], y[n] -> out[n]
double sum = 0.0;
for (unsigned n = 0; n < signal_len; n++) {
    sum += signal[n];
    if (n >= window_len)
        sum -= signal[n - window_len];
    out[n] = sum / (double)window_len;
}
```

**API:**

```c
void MD_moving_average(const double *signal, unsigned signal_len,
                       unsigned window_len, double *out);
```

**Visuals** — response and magnitude spectrum:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="moving_average_response.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="moving_average_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

**Quick example:**

\snippet fir_convolution.c moving-average

---

## General FIR filter

With arbitrary coefficients \f$b[k]\f$ (length \f$M\f$), causal FIR filtering is:

\f[
y[n] = \sum_{k=0}^{M-1} b[k]\,x[n-k], \quad n = 0,1,\ldots,N-1
\f]

This is the direct form used for most textbook FIR designs.

**Reading the formula in C:**

```c
// b[k] -> coeffs[k], x[n-k] -> signal[n-k], y[n] -> out[n]
for (unsigned n = 0; n < signal_len; n++) {
    double acc = 0.0;
    for (unsigned k = 0; k < num_taps; k++) {
        if (k > n) break;    // zero-padded startup
        acc += coeffs[k] * signal[n - k];
    }
    out[n] = acc;
}
```

**API:**

```c
void MD_fir_filter(const double *signal, unsigned signal_len,
                   const double *coeffs, unsigned num_taps,
                   double *out);
```

**Visuals** — response and magnitude spectrum:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="fir_general_response.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="fir_general_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

**Quick example:**

\snippet fir_convolution.c fir-filter

---

## FFT overlap-add fast convolution

For long filters, direct convolution is \f$O(NM)\f$. [Overlap-add](https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method) computes the
same full linear convolution in blocks using FFTs:

\f[
Y_b(k) = X_b(k)\,H(k), \quad
y_b[n] = \mathrm{IFFT}(Y_b(k))
\f]

Each block's time-domain result is added into the output with overlap
between adjacent blocks, after an [inverse FFT (IFFT)](https://en.wikipedia.org/wiki/Fast_Fourier_transform#Inverse_transform) per block.

**Reading the algorithm in C:**

```c
// 1) Zero-pad one input block to nfft and FFT it
// 2) Multiply by precomputed FFT(kernel)
// 3) IFFT back to time domain
// 4) Add (overlap-add) into full output buffer
for (unsigned start = 0; start < signal_len; start += block_len) {
    // load block, FFT, multiply in frequency domain, IFFT, overlap-add
}
```

**API:**

```c
void MD_convolution_fft_ola(const double *signal, unsigned signal_len,
                            const double *kernel, unsigned kernel_len,
                            double *out);
```

**Visuals** — response and magnitude spectrum:

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="conv_fft_ola_response.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="conv_fft_ola_spectrum.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

**Quick example:**

\snippet fir_convolution.c fft-convolution

---

## Quick comparison

| Method | Output length | Typical cost | Best use |
|--------|---------------|--------------|----------|
| Time-domain convolution | \f$N+M-1\f$ | \f$O(NM)\f$ | Teaching and short kernels |
| Moving average | \f$N\f$ | \f$O(N)\f$ (running sum) | First FIR low-pass example |
| General FIR | \f$N\f$ | \f$O(NM)\f$ | Arbitrary FIR taps |
| FFT overlap-add | \f$N+M-1\f$ | \f$O(N\log N)\f$ blocks | Long kernels / fast offline convolution |

Use direct methods to build intuition, then compare with overlap-add to see
why FFT-based convolution matters for longer filters.

## Further reading

- [Convolution](https://en.wikipedia.org/wiki/Convolution) -- Wikipedia
- [Finite impulse response](https://en.wikipedia.org/wiki/Finite_impulse_response) -- Wikipedia
- [Overlap-add method](https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method) -- Wikipedia

## API reference

- MD_convolution_num_samples() -- compute full-convolution output length
- MD_convolution_time() -- direct time-domain full convolution
- MD_moving_average() -- causal moving-average FIR filter
- MD_fir_filter() -- causal arbitrary-coefficient FIR filter
- MD_convolution_fft_ola() -- full convolution via FFT overlap-add
