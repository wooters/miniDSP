/**
 * @file minidsp_fir.c
 * @brief FIR filtering and convolution: direct time-domain methods and
 *        FFT overlap-add fast convolution.
 */

#include "minidsp.h"

/** Return the next power-of-two >= n (n must be > 0). */
static unsigned next_pow2(unsigned n)
{
    unsigned p = 1;
    while (p < n) p <<= 1;
    return p;
}

unsigned MD_convolution_num_samples(unsigned signal_len, unsigned kernel_len)
{
    assert(signal_len > 0);
    assert(kernel_len > 0);
    return signal_len + kernel_len - 1;
}

void MD_convolution_time(const double *signal, unsigned signal_len,
                         const double *kernel, unsigned kernel_len,
                         double *out)
{
    assert(signal != nullptr);
    assert(kernel != nullptr);
    assert(out != nullptr);
    assert(signal_len > 0);
    assert(kernel_len > 0);

    unsigned out_len = MD_convolution_num_samples(signal_len, kernel_len);
    for (unsigned n = 0; n < out_len; n++) {
        double acc = 0.0;
        for (unsigned k = 0; k < kernel_len; k++) {
            if (k > n) break;                /* signal index would be negative */
            unsigned si = n - k;
            if (si >= signal_len) continue;  /* past end of input signal */
            acc += signal[si] * kernel[k];
        }
        out[n] = acc;
    }
}

void MD_moving_average(const double *signal, unsigned signal_len,
                       unsigned window_len, double *out)
{
    assert(signal != nullptr);
    assert(out != nullptr);
    assert(signal_len > 0);
    assert(window_len > 0);

    double sum = 0.0;
    double inv_win = 1.0 / (double)window_len;

    for (unsigned n = 0; n < signal_len; n++) {
        sum += signal[n];
        if (n >= window_len) {
            sum -= signal[n - window_len];
        }
        out[n] = sum * inv_win;
    }
}

void MD_fir_filter(const double *signal, unsigned signal_len,
                   const double *coeffs, unsigned num_taps,
                   double *out)
{
    assert(signal != nullptr);
    assert(coeffs != nullptr);
    assert(out != nullptr);
    assert(signal_len > 0);
    assert(num_taps > 0);

    for (unsigned n = 0; n < signal_len; n++) {
        double acc = 0.0;
        for (unsigned k = 0; k < num_taps; k++) {
            if (k > n) break;  /* zero-padded causal startup */
            acc += coeffs[k] * signal[n - k];
        }
        out[n] = acc;
    }
}

void MD_convolution_fft_ola(const double *signal, unsigned signal_len,
                            const double *kernel, unsigned kernel_len,
                            double *out)
{
    assert(signal != nullptr);
    assert(kernel != nullptr);
    assert(out != nullptr);
    assert(signal_len > 0);
    assert(kernel_len > 0);

    unsigned out_len = MD_convolution_num_samples(signal_len, kernel_len);
    memset(out, 0, out_len * sizeof(double));

    unsigned nfft = next_pow2(2U * kernel_len);
    if (nfft < 64U) nfft = 64U;
    unsigned block_len = nfft - kernel_len + 1U;

    unsigned num_bins = nfft / 2U + 1U;

    double *h_time = calloc(nfft, sizeof(double));
    double *x_time = calloc(nfft, sizeof(double));
    double *y_time = calloc(nfft, sizeof(double));
    fftw_complex *H = fftw_alloc_complex(num_bins);
    fftw_complex *X = fftw_alloc_complex(num_bins);
    fftw_complex *Y = fftw_alloc_complex(num_bins);
    assert(h_time != nullptr);
    assert(x_time != nullptr);
    assert(y_time != nullptr);
    assert(H != nullptr);
    assert(X != nullptr);
    assert(Y != nullptr);

    memcpy(h_time, kernel, kernel_len * sizeof(double));

    fftw_plan plan_h = fftw_plan_dft_r2c_1d((int)nfft, h_time, H, FFTW_ESTIMATE);
    fftw_plan plan_x = fftw_plan_dft_r2c_1d((int)nfft, x_time, X, FFTW_ESTIMATE);
    fftw_plan plan_y = fftw_plan_dft_c2r_1d((int)nfft, Y, y_time, FFTW_ESTIMATE);
    assert(plan_h != nullptr);
    assert(plan_x != nullptr);
    assert(plan_y != nullptr);

    fftw_execute(plan_h);  /* kernel FFT once */

    for (unsigned start = 0; start < signal_len; start += block_len) {
        unsigned in_len = signal_len - start;
        if (in_len > block_len) in_len = block_len;

        memset(x_time, 0, nfft * sizeof(double));
        memcpy(x_time, signal + start, in_len * sizeof(double));

        fftw_execute(plan_x);

        for (unsigned k = 0; k < num_bins; k++) {
            double xr = creal(X[k]);
            double xi = cimag(X[k]);
            double hr = creal(H[k]);
            double hi = cimag(H[k]);
            Y[k] = (xr * hr - xi * hi) + I * (xr * hi + xi * hr);
        }

        fftw_execute(plan_y);

        double norm = 1.0 / (double)nfft;
        for (unsigned n = 0; n < nfft; n++) {
            unsigned oi = start + n;
            if (oi >= out_len) break;
            out[oi] += y_time[n] * norm;
        }
    }

    fftw_destroy_plan(plan_h);
    fftw_destroy_plan(plan_x);
    fftw_destroy_plan(plan_y);
    fftw_free(Y);
    fftw_free(X);
    fftw_free(H);
    free(y_time);
    free(x_time);
    free(h_time);
}
