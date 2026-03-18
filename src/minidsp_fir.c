/**
 * @file minidsp_fir.c
 * @brief FIR filtering and convolution: direct time-domain methods and
 *        FFT overlap-add fast convolution.
 */

#include "minidsp.h"
#include "minidsp_internal.h"

/** Return the next power-of-two >= n (n must be > 0). */
static unsigned next_pow2(unsigned n)
{
    unsigned p = 1;
    while (p < n) p <<= 1;
    return p;
}

unsigned MD_convolution_num_samples(unsigned signal_len, unsigned kernel_len)
{
    MD_CHECK(signal_len > 0, MD_ERR_INVALID_SIZE, "signal_len must be > 0", 0);
    MD_CHECK(kernel_len > 0, MD_ERR_INVALID_SIZE, "kernel_len must be > 0", 0);
    return signal_len + kernel_len - 1;
}

void MD_convolution_time(const double *signal, unsigned signal_len,
                         const double *kernel, unsigned kernel_len,
                         double *out)
{
    MD_CHECK_VOID(signal != NULL, MD_ERR_NULL_POINTER, "signal must not be NULL");
    MD_CHECK_VOID(kernel != NULL, MD_ERR_NULL_POINTER, "kernel must not be NULL");
    MD_CHECK_VOID(out != NULL, MD_ERR_NULL_POINTER, "out must not be NULL");
    MD_CHECK_VOID(signal_len > 0, MD_ERR_INVALID_SIZE, "signal_len must be > 0");
    MD_CHECK_VOID(kernel_len > 0, MD_ERR_INVALID_SIZE, "kernel_len must be > 0");

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
    MD_CHECK_VOID(signal != NULL, MD_ERR_NULL_POINTER, "signal must not be NULL");
    MD_CHECK_VOID(out != NULL, MD_ERR_NULL_POINTER, "out must not be NULL");
    MD_CHECK_VOID(signal_len > 0, MD_ERR_INVALID_SIZE, "signal_len must be > 0");
    MD_CHECK_VOID(window_len > 0, MD_ERR_INVALID_SIZE, "window_len must be > 0");

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
    MD_CHECK_VOID(signal != NULL, MD_ERR_NULL_POINTER, "signal must not be NULL");
    MD_CHECK_VOID(coeffs != NULL, MD_ERR_NULL_POINTER, "coeffs must not be NULL");
    MD_CHECK_VOID(out != NULL, MD_ERR_NULL_POINTER, "out must not be NULL");
    MD_CHECK_VOID(signal_len > 0, MD_ERR_INVALID_SIZE, "signal_len must be > 0");
    MD_CHECK_VOID(num_taps > 0, MD_ERR_INVALID_SIZE, "num_taps must be > 0");

    for (unsigned n = 0; n < signal_len; n++) {
        double acc = 0.0;
        for (unsigned k = 0; k < num_taps; k++) {
            if (k > n) break;  /* zero-padded causal startup */
            acc += coeffs[k] * signal[n - k];
        }
        out[n] = acc;
    }
}

void MD_design_lowpass_fir(double *coeffs, unsigned num_taps,
                           double cutoff_freq, double sample_rate,
                           double kaiser_beta)
{
    MD_CHECK_VOID(coeffs != NULL, MD_ERR_NULL_POINTER, "coeffs must not be NULL");
    MD_CHECK_VOID(num_taps > 0, MD_ERR_INVALID_SIZE, "num_taps must be > 0");
    MD_CHECK_VOID(sample_rate > 0.0, MD_ERR_INVALID_RANGE, "sample_rate must be > 0");
    MD_CHECK_VOID(cutoff_freq > 0.0, MD_ERR_INVALID_RANGE, "cutoff_freq must be > 0");
    MD_CHECK_VOID(cutoff_freq < sample_rate / 2.0, MD_ERR_INVALID_RANGE, "cutoff_freq must be < sample_rate/2");

    double fc = cutoff_freq / sample_rate;  /* normalized cutoff */
    double center = (double)(num_taps - 1) / 2.0;

    /* Generate Kaiser window */
    double *kaiser = malloc(num_taps * sizeof(double));
    MD_CHECK_VOID(kaiser != NULL, MD_ERR_ALLOC_FAILED, "kaiser window allocation failed");
    MD_Gen_Kaiser_Win(kaiser, num_taps, kaiser_beta);

    /* Windowed sinc */
    double sum = 0.0;
    for (unsigned i = 0; i < num_taps; i++) {
        double x = (double)i - center;
        coeffs[i] = 2.0 * fc * MD_sinc(2.0 * fc * x) * kaiser[i];
        sum += coeffs[i];
    }

    free(kaiser);

    /* Normalize for unity DC gain */
    if (fabs(sum) > 1e-30) {
        for (unsigned i = 0; i < num_taps; i++) {
            coeffs[i] /= sum;
        }
    }
}

void MD_convolution_fft_ola(const double *signal, unsigned signal_len,
                            const double *kernel, unsigned kernel_len,
                            double *out)
{
    MD_CHECK_VOID(signal != NULL, MD_ERR_NULL_POINTER, "signal must not be NULL");
    MD_CHECK_VOID(kernel != NULL, MD_ERR_NULL_POINTER, "kernel must not be NULL");
    MD_CHECK_VOID(out != NULL, MD_ERR_NULL_POINTER, "out must not be NULL");
    MD_CHECK_VOID(signal_len > 0, MD_ERR_INVALID_SIZE, "signal_len must be > 0");
    MD_CHECK_VOID(kernel_len > 0, MD_ERR_INVALID_SIZE, "kernel_len must be > 0");

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
    MD_CHECK_VOID(h_time != NULL, MD_ERR_ALLOC_FAILED, "h_time allocation failed");
    MD_CHECK_VOID(x_time != NULL, MD_ERR_ALLOC_FAILED, "x_time allocation failed");
    MD_CHECK_VOID(y_time != NULL, MD_ERR_ALLOC_FAILED, "y_time allocation failed");
    MD_CHECK_VOID(H != NULL, MD_ERR_ALLOC_FAILED, "H allocation failed");
    MD_CHECK_VOID(X != NULL, MD_ERR_ALLOC_FAILED, "X allocation failed");
    MD_CHECK_VOID(Y != NULL, MD_ERR_ALLOC_FAILED, "Y allocation failed");

    memcpy(h_time, kernel, kernel_len * sizeof(double));

    fftw_plan plan_h = fftw_plan_dft_r2c_1d((int)nfft, h_time, H, FFTW_ESTIMATE);
    fftw_plan plan_x = fftw_plan_dft_r2c_1d((int)nfft, x_time, X, FFTW_ESTIMATE);
    fftw_plan plan_y = fftw_plan_dft_c2r_1d((int)nfft, Y, y_time, FFTW_ESTIMATE);
    MD_CHECK_VOID(plan_h != NULL, MD_ERR_ALLOC_FAILED, "plan_h allocation failed");
    MD_CHECK_VOID(plan_x != NULL, MD_ERR_ALLOC_FAILED, "plan_x allocation failed");
    MD_CHECK_VOID(plan_y != NULL, MD_ERR_ALLOC_FAILED, "plan_y allocation failed");

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
