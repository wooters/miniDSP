/**
 * @file fir_convolution.c
 * @brief Example: FIR filtering and convolution in time and frequency domains.
 *
 * Demonstrates:
 *   1. MD_convolution_time() -- direct full convolution
 *   2. MD_moving_average() -- causal moving-average filter
 *   3. MD_fir_filter() -- causal arbitrary-coefficient FIR
 *   4. MD_convolution_fft_ola() -- fast full convolution (overlap-add FFT)
 *
 * Build and run (from repo root):
 *   make -C examples fir_convolution
 *   ./examples/fir_convolution
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"

static void print_first(const char *name, const double *x, unsigned n, unsigned max_show)
{
    unsigned show = (n < max_show) ? n : max_show;
    printf("%s (first %u):", name, show);
    for (unsigned i = 0; i < show; i++) {
        printf(" % .6f", x[i]);
    }
    printf("\n");
}

int main(void)
{
    const unsigned N = 1024;
    const double sample_rate = 16000.0;

    double *signal = malloc(N * sizeof(double));
    if (!signal) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    /* Build a test input: two tones plus a small DC offset. */
    for (unsigned i = 0; i < N; i++) {
        signal[i] = 0.8 * sin(2.0 * M_PI * 440.0 * (double)i / sample_rate)
                  + 0.2 * sin(2.0 * M_PI * 2000.0 * (double)i / sample_rate)
                  + 0.05;
    }

    /* ---------------------------------------------------------------
     * 1) Time-domain full convolution
     * ---------------------------------------------------------------*/
    double kernel[] = {0.2, 0.6, 0.2};
    unsigned kernel_len = sizeof(kernel) / sizeof(kernel[0]);
    unsigned conv_len = MD_convolution_num_samples(N, kernel_len);
    double *y_conv_time = malloc(conv_len * sizeof(double));
    double *y_conv_fft = malloc(conv_len * sizeof(double));
    if (!y_conv_time || !y_conv_fft) {
        fprintf(stderr, "allocation failed\n");
        free(y_conv_fft);
        free(y_conv_time);
        free(signal);
        return 1;
    }

    //! [convolution-time]
    MD_convolution_time(signal, N, kernel, kernel_len, y_conv_time);
    //! [convolution-time]

    /* ---------------------------------------------------------------
     * 2) Moving-average filter
     * ---------------------------------------------------------------*/
    double *y_ma = malloc(N * sizeof(double));
    if (!y_ma) {
        fprintf(stderr, "allocation failed\n");
        free(y_conv_fft);
        free(y_conv_time);
        free(signal);
        return 1;
    }

    //! [moving-average]
    MD_moving_average(signal, N, 8, y_ma);
    //! [moving-average]

    /* ---------------------------------------------------------------
     * 3) General FIR filter
     * ---------------------------------------------------------------*/
    double coeffs[] = {0.05, 0.10, 0.20, 0.30, 0.20, 0.10, 0.05};
    unsigned num_taps = sizeof(coeffs) / sizeof(coeffs[0]);
    double *y_fir = malloc(N * sizeof(double));
    if (!y_fir) {
        fprintf(stderr, "allocation failed\n");
        free(y_ma);
        free(y_conv_fft);
        free(y_conv_time);
        free(signal);
        return 1;
    }

    //! [fir-filter]
    MD_fir_filter(signal, N, coeffs, num_taps, y_fir);
    //! [fir-filter]

    /* ---------------------------------------------------------------
     * 4) Fast convolution via overlap-add FFT
     * ---------------------------------------------------------------*/
    //! [fft-convolution]
    MD_convolution_fft_ola(signal, N, kernel, kernel_len, y_conv_fft);
    //! [fft-convolution]

    /* Compare direct and FFT full convolution results. */
    double max_err = 0.0;
    for (unsigned i = 0; i < conv_len; i++) {
        double e = fabs(y_conv_time[i] - y_conv_fft[i]);
        if (e > max_err) max_err = e;
    }

    print_first("Input", signal, N, 8);
    print_first("Time convolution", y_conv_time, conv_len, 8);
    print_first("Moving average", y_ma, N, 8);
    print_first("General FIR", y_fir, N, 8);
    print_first("FFT OLA convolution", y_conv_fft, conv_len, 8);
    printf("Max |time_conv - fft_conv| = %.3e\n", max_err);

    MD_shutdown();
    free(y_fir);
    free(y_ma);
    free(y_conv_fft);
    free(y_conv_time);
    free(signal);
    return 0;
}
