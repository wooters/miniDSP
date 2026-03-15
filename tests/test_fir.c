/**
 * @file test_fir.c
 * @brief Tests for FIR filters, convolution, Bessel, sinc, Kaiser, and
 *        lowpass FIR design.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * Tests for FIR filtering / convolution
 * -----------------------------------------------------------------------*/

static int test_convolution_num_samples(void)
{
    int ok = 1;
    ok &= (MD_convolution_num_samples(1, 1) == 1);
    ok &= (MD_convolution_num_samples(4, 3) == 6);
    ok &= (MD_convolution_num_samples(256, 17) == 272);
    return ok;
}

/** Small hand-computed full convolution. */
static int test_convolution_time_known_small(void)
{
    double x[] = {1.0, 2.0, 3.0};
    double h[] = {1.0, 1.0};
    double y[4];
    double expected[] = {1.0, 3.0, 5.0, 3.0};

    MD_convolution_time(x, 3, h, 2, y);

    int ok = 1;
    for (unsigned i = 0; i < 4; i++) {
        ok &= approx_equal(y[i], expected[i], 1e-15);
    }
    return ok;
}

/** Convolution with a unit impulse should reproduce the input. */
static int test_convolution_time_impulse_identity(void)
{
    double x[] = {-1.0, 0.5, 2.0, -3.0, 4.5};
    double h[] = {1.0};
    double y[5];

    MD_convolution_time(x, 5, h, 1, y);

    int ok = 1;
    for (unsigned i = 0; i < 5; i++) {
        ok &= approx_equal(y[i], x[i], 1e-15);
    }
    return ok;
}

/** Known 2-tap FIR coefficients should match hand-computed values. */
static int test_fir_filter_known_taps(void)
{
    double x[] = {1.0, 2.0, 3.0, 4.0};
    double b[] = {0.5, 0.5};
    double y[4];
    double expected[] = {0.5, 1.5, 2.5, 3.5};

    MD_fir_filter(x, 4, b, 2, y);

    int ok = 1;
    for (unsigned i = 0; i < 4; i++) {
        ok &= approx_equal(y[i], expected[i], 1e-15);
    }
    return ok;
}

/** FIR filter with single coefficient scales the signal. */
static int test_fir_filter_single_tap(void)
{
    double x[] = {1.0, 2.0, 3.0, 4.0};
    double b[] = {0.5};
    double y[4];
    MD_fir_filter(x, 4, b, 1, y);
    int ok = 1;
    for (unsigned i = 0; i < 4; i++) {
        ok &= approx_equal(y[i], x[i] * 0.5, 1e-15);
    }
    return ok;
}

/** Step response of a causal moving average with fixed denominator. */
static int test_moving_average_step_response(void)
{
    double x[] = {1, 1, 1, 1, 1, 1, 1, 1};
    double y[8];
    double expected[] = {0.25, 0.50, 0.75, 1.0, 1.0, 1.0, 1.0, 1.0};

    MD_moving_average(x, 8, 4, y);

    int ok = 1;
    for (unsigned i = 0; i < 8; i++) {
        ok &= approx_equal(y[i], expected[i], 1e-15);
    }
    return ok;
}

/** Moving average should match FIR with boxcar coefficients (1/M each). */
static int test_moving_average_matches_boxcar_fir(void)
{
    unsigned N = 256;
    unsigned M = 9;

    double *x = malloc(N * sizeof(double));
    double *y_ma = malloc(N * sizeof(double));
    double *y_fir = malloc(N * sizeof(double));
    double *b = malloc(M * sizeof(double));
    if (!x || !y_ma || !y_fir || !b) {
        free(b);
        free(y_fir);
        free(y_ma);
        free(x);
        return 0;
    }

    for (unsigned i = 0; i < N; i++) {
        x[i] = sin(2.0 * M_PI * 11.0 * (double)i / (double)N)
             + 0.3 * cos(2.0 * M_PI * 37.0 * (double)i / (double)N);
    }
    for (unsigned k = 0; k < M; k++) {
        b[k] = 1.0 / (double)M;
    }

    MD_moving_average(x, N, M, y_ma);
    MD_fir_filter(x, N, b, M, y_fir);

    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        ok &= approx_equal(y_ma[i], y_fir[i], 1e-12);
    }

    free(b);
    free(y_fir);
    free(y_ma);
    free(x);
    return ok;
}

/** Moving average with width=1 is identity (output equals input). */
static int test_moving_average_width_one(void)
{
    double x[] = {1.0, -2.0, 3.0, -4.0, 5.0};
    double y[5];
    MD_moving_average(x, 5, 1, y);
    int ok = 1;
    for (unsigned i = 0; i < 5; i++) {
        ok &= approx_equal(y[i], x[i], 1e-15);
    }
    return ok;
}

/** FFT overlap-add should match direct time-domain convolution. */
static int test_convolution_fft_ola_matches_time(void)
{
    unsigned N = 2048;
    unsigned M = 257;
    unsigned out_len = MD_convolution_num_samples(N, M);

    double *x = malloc(N * sizeof(double));
    double *h = malloc(M * sizeof(double));
    double *y_time = malloc(out_len * sizeof(double));
    double *y_fft = malloc(out_len * sizeof(double));
    if (!x || !h || !y_time || !y_fft) {
        free(y_fft);
        free(y_time);
        free(h);
        free(x);
        return 0;
    }

    for (unsigned i = 0; i < N; i++) {
        x[i] = 0.6 * sin(2.0 * M_PI * 37.0 * (double)i / (double)N)
             + 0.2 * cos(2.0 * M_PI * 113.0 * (double)i / (double)N);
    }
    for (unsigned k = 0; k < M; k++) {
        /* Finite decaying kernel */
        h[k] = exp(-0.01 * (double)k) * cos(2.0 * M_PI * 0.02 * (double)k);
    }

    MD_convolution_time(x, N, h, M, y_time);
    MD_convolution_fft_ola(x, N, h, M, y_fft);

    double max_err = 0.0;
    for (unsigned i = 0; i < out_len; i++) {
        double err = fabs(y_time[i] - y_fft[i]);
        if (err > max_err) max_err = err;
    }

    free(y_fft);
    free(y_time);
    free(h);
    free(x);
    return (max_err <= 1e-8);
}

/** FFT overlap-add should remain accurate across length combinations. */
static int test_convolution_fft_ola_different_lengths(void)
{
    const unsigned signal_lens[] = {63, 500, 409};
    const unsigned kernel_lens[] = {5, 80, 257};
    const unsigned num_cases = 3;

    int ok = 1;
    for (unsigned c = 0; c < num_cases; c++) {
        unsigned N = signal_lens[c];
        unsigned M = kernel_lens[c];
        unsigned out_len = MD_convolution_num_samples(N, M);

        double *x = malloc(N * sizeof(double));
        double *h = malloc(M * sizeof(double));
        double *y_time = malloc(out_len * sizeof(double));
        double *y_fft = malloc(out_len * sizeof(double));
        if (!x || !h || !y_time || !y_fft) {
            free(y_fft);
            free(y_time);
            free(h);
            free(x);
            return 0;
        }

        for (unsigned i = 0; i < N; i++) {
            x[i] = sin(2.0 * M_PI * (3.0 + (double)c) * (double)i / (double)N)
                 + 0.1 * (double)((int)(i % 7) - 3);
        }
        for (unsigned k = 0; k < M; k++) {
            h[k] = 1.0 / (double)(k + 1);
        }

        MD_convolution_time(x, N, h, M, y_time);
        MD_convolution_fft_ola(x, N, h, M, y_fft);

        for (unsigned i = 0; i < out_len; i++) {
            if (!approx_equal(y_time[i], y_fft[i], 1e-8)) {
                ok = 0;
                break;
            }
        }

        free(y_fft);
        free(y_time);
        free(h);
        free(x);
        if (!ok) break;
    }

    return ok;
}

/** FFT overlap-add with length-1 unit kernel reproduces input. */
static int test_convolution_fft_ola_impulse_identity(void)
{
    double x[] = {1.0, -2.0, 3.0, -4.0, 5.0, -6.0, 7.0, -8.0};
    double h[] = {1.0};
    unsigned out_len = MD_convolution_num_samples(8, 1);
    double y[8];
    MD_convolution_fft_ola(x, 8, h, 1, y);
    int ok = 1;
    for (unsigned i = 0; i < 8; i++) {
        ok &= approx_equal(y[i], x[i], 1e-10);
    }
    (void)out_len;
    return ok;
}

/* -----------------------------------------------------------------------
 * MD_bessel_i0 tests
 * -----------------------------------------------------------------------*/

static int test_bessel_i0_at_zero(void)
{
    return approx_equal(MD_bessel_i0(0.0), 1.0, 1e-15);
}

static int test_bessel_i0_known_x1(void)
{
    return approx_equal(MD_bessel_i0(1.0), 1.2660658777520084, 1e-10);
}

static int test_bessel_i0_known_x5(void)
{
    return approx_equal(MD_bessel_i0(5.0), 27.239871823604442, 1e-6);
}

static int test_bessel_i0_monotonic(void)
{
    double prev = MD_bessel_i0(0.0);
    double xs[] = {1.0, 2.0, 3.0, 4.0, 5.0};
    for (unsigned i = 0; i < 5; i++) {
        double cur = MD_bessel_i0(xs[i]);
        if (cur <= prev) return 0;
        prev = cur;
    }
    return 1;
}

/** I0 is an even function: I0(-x) == I0(x). */
static int test_bessel_i0_negative_symmetry(void)
{
    double xs[] = {1.0, 2.5, 5.0, 10.0};
    for (unsigned i = 0; i < 4; i++) {
        if (!approx_equal(MD_bessel_i0(-xs[i]), MD_bessel_i0(xs[i]), 1e-12))
            return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * MD_sinc tests
 * -----------------------------------------------------------------------*/

static int test_sinc_at_zero(void)
{
    return approx_equal(MD_sinc(0.0), 1.0, 1e-15);
}

static int test_sinc_integer_zeros(void)
{
    int ints[] = {-3, -2, -1, 1, 2, 3};
    for (unsigned i = 0; i < 6; i++) {
        if (!approx_equal(MD_sinc((double)ints[i]), 0.0, 1e-15))
            return 0;
    }
    return 1;
}

static int test_sinc_half(void)
{
    return approx_equal(MD_sinc(0.5), 2.0 / M_PI, 1e-12);
}

static int test_sinc_near_zero_threshold(void)
{
    return approx_equal(MD_sinc(1e-13), 1.0, 1e-15);
}

/** Sinc is an even function: sinc(-x) == sinc(x). */
static int test_sinc_symmetry(void)
{
    double xs[] = {0.5, 1.5, 2.7, 0.01, 3.14};
    for (unsigned i = 0; i < 5; i++) {
        if (!approx_equal(MD_sinc(-xs[i]), MD_sinc(xs[i]), 1e-15))
            return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * MD_Gen_Kaiser_Win tests
 * -----------------------------------------------------------------------*/

static int test_kaiser_single_sample(void)
{
    double w;
    MD_Gen_Kaiser_Win(&w, 1, 10.0);
    return approx_equal(w, 1.0, 1e-15);
}

static int test_kaiser_symmetry(void)
{
    unsigned n = 128;
    double w[128];
    MD_Gen_Kaiser_Win(w, n, 10.0);
    for (unsigned i = 0; i < n / 2; i++) {
        if (!approx_equal(w[i], w[n - 1 - i], 1e-12))
            return 0;
    }
    return 1;
}

static int test_kaiser_tapered_ends(void)
{
    unsigned n = 128;
    double w[128];
    MD_Gen_Kaiser_Win(w, n, 10.0);
    /* Ends should be positive but less than center */
    if (w[0] <= 0.0) return 0;
    if (w[0] >= w[n / 2]) return 0;
    return 1;
}

static int test_kaiser_peak_at_center(void)
{
    unsigned n = 129;
    double w[129];
    MD_Gen_Kaiser_Win(w, n, 10.0);
    unsigned center = n / 2;
    /* Center should be the maximum */
    for (unsigned i = 0; i < n; i++) {
        if (w[i] > w[center] + 1e-12) return 0;
    }
    return 1;
}

static int test_kaiser_beta_comparison(void)
{
    unsigned n = 128;
    double w_low[128], w_high[128];
    MD_Gen_Kaiser_Win(w_low, n, 5.0);
    MD_Gen_Kaiser_Win(w_high, n, 14.0);
    /* Higher beta → more tapered (smaller end/center ratio) */
    double ratio_low = w_low[0] / w_low[n / 2];
    double ratio_high = w_high[0] / w_high[n / 2];
    return ratio_high < ratio_low;
}

/** Kaiser window with beta=0 should equal a rectangular window (all ones). */
static int test_kaiser_beta_zero(void)
{
    unsigned n = 64;
    double w[64];
    MD_Gen_Kaiser_Win(w, n, 0.0);
    for (unsigned i = 0; i < n; i++) {
        if (!approx_equal(w[i], 1.0, 1e-12)) return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * MD_design_lowpass_fir tests
 * -----------------------------------------------------------------------*/

static int test_lowpass_fir_dc_gain(void)
{
    unsigned taps = 65;
    double h[65];
    MD_design_lowpass_fir(h, taps, 4000.0, 48000.0, 10.0);
    double sum = 0.0;
    for (unsigned i = 0; i < taps; i++) sum += h[i];
    return approx_equal(sum, 1.0, 1e-12);
}

static int test_lowpass_fir_symmetry(void)
{
    unsigned taps = 65;
    double h[65];
    MD_design_lowpass_fir(h, taps, 4000.0, 48000.0, 10.0);
    for (unsigned i = 0; i < taps / 2; i++) {
        if (!approx_equal(h[i], h[taps - 1 - i], 1e-12))
            return 0;
    }
    return 1;
}

static int test_lowpass_fir_passband(void)
{
    /* A 1 kHz tone through a 4 kHz LPF should pass with < 0.1 dB loss */
    unsigned sr = 48000;
    unsigned N = 4096;
    double *sig = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));
    unsigned taps = 65;
    double h[65];

    MD_sine_wave(sig, N, 1.0, 1000.0, (double)sr);
    MD_design_lowpass_fir(h, taps, 4000.0, (double)sr, 10.0);
    MD_fir_filter(sig, N, h, taps, out);

    /* Measure RMS of steady-state output (skip startup transient) */
    double rms_in = MD_rms(sig + taps, N - taps);
    double rms_out = MD_rms(out + taps, N - taps);
    double db_diff = 20.0 * log10(rms_out / rms_in);

    free(sig);
    free(out);
    return fabs(db_diff) < 0.1;
}

static int test_lowpass_fir_stopband(void)
{
    /* A 10 kHz tone through a 4 kHz LPF should be attenuated > 60 dB */
    unsigned sr = 48000;
    unsigned N = 4096;
    double *sig = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));
    unsigned taps = 65;
    double h[65];

    MD_sine_wave(sig, N, 1.0, 10000.0, (double)sr);
    MD_design_lowpass_fir(h, taps, 4000.0, (double)sr, 10.0);
    MD_fir_filter(sig, N, h, taps, out);

    double rms_in = MD_rms(sig + taps, N - taps);
    double rms_out = MD_rms(out + taps, N - taps);
    double attenuation_db = 20.0 * log10(rms_out / rms_in);

    free(sig);
    free(out);
    return attenuation_db < -60.0;
}

/** Near-Nyquist cutoff should pass a high-frequency tone with minimal attenuation. */
static int test_lowpass_fir_near_nyquist(void)
{
    unsigned sr = 48000;
    unsigned N = 4096;
    double *sig = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));
    unsigned taps = 65;
    double h[65];

    /* 10 kHz tone through a 20 kHz LPF (Nyquist at 24 kHz) */
    MD_sine_wave(sig, N, 1.0, 10000.0, (double)sr);
    MD_design_lowpass_fir(h, taps, 20000.0, (double)sr, 10.0);
    MD_fir_filter(sig, N, h, taps, out);

    double rms_in = MD_rms(sig + taps, N - taps);
    double rms_out = MD_rms(out + taps, N - taps);
    double db_diff = 20.0 * log10(rms_out / rms_in);

    free(sig);
    free(out);
    return fabs(db_diff) < 0.5;
}

/** Very low cutoff should heavily attenuate a mid-range tone. */
static int test_lowpass_fir_very_low_cutoff(void)
{
    unsigned sr = 44100;
    unsigned N = 8192;
    unsigned taps = 257;
    double *sig = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));
    double *h   = malloc(taps * sizeof(double));

    /* 5 kHz tone through a 100 Hz LPF with enough taps for sharp rolloff */
    MD_sine_wave(sig, N, 1.0, 5000.0, (double)sr);
    MD_design_lowpass_fir(h, taps, 100.0, (double)sr, 10.0);
    MD_fir_filter(sig, N, h, taps, out);

    double rms_in = MD_rms(sig + taps, N - taps);
    double rms_out = MD_rms(out + taps, N - taps);
    double attenuation_db = 20.0 * log10(rms_out / rms_in);

    free(h);
    free(sig);
    free(out);
    return attenuation_db < -30.0;
}

/* -----------------------------------------------------------------------
 * Runner — called from main() in test_minidsp.c
 * -----------------------------------------------------------------------*/

void run_fir_tests(void)
{
    printf("\n--- FIR filters / convolution ---\n");
    RUN_TEST(test_convolution_num_samples);
    RUN_TEST(test_convolution_time_known_small);
    RUN_TEST(test_convolution_time_impulse_identity);
    RUN_TEST(test_fir_filter_known_taps);
    RUN_TEST(test_fir_filter_single_tap);
    RUN_TEST(test_moving_average_step_response);
    RUN_TEST(test_moving_average_matches_boxcar_fir);
    RUN_TEST(test_moving_average_width_one);
    RUN_TEST(test_convolution_fft_ola_matches_time);
    RUN_TEST(test_convolution_fft_ola_different_lengths);
    RUN_TEST(test_convolution_fft_ola_impulse_identity);

    printf("\n--- MD_bessel_i0 ---\n");
    RUN_TEST(test_bessel_i0_at_zero);
    RUN_TEST(test_bessel_i0_known_x1);
    RUN_TEST(test_bessel_i0_known_x5);
    RUN_TEST(test_bessel_i0_monotonic);
    RUN_TEST(test_bessel_i0_negative_symmetry);

    printf("\n--- MD_sinc ---\n");
    RUN_TEST(test_sinc_at_zero);
    RUN_TEST(test_sinc_integer_zeros);
    RUN_TEST(test_sinc_half);
    RUN_TEST(test_sinc_near_zero_threshold);
    RUN_TEST(test_sinc_symmetry);

    printf("\n--- MD_Gen_Kaiser_Win ---\n");
    RUN_TEST(test_kaiser_single_sample);
    RUN_TEST(test_kaiser_symmetry);
    RUN_TEST(test_kaiser_tapered_ends);
    RUN_TEST(test_kaiser_peak_at_center);
    RUN_TEST(test_kaiser_beta_comparison);
    RUN_TEST(test_kaiser_beta_zero);

    printf("\n--- MD_design_lowpass_fir ---\n");
    RUN_TEST(test_lowpass_fir_dc_gain);
    RUN_TEST(test_lowpass_fir_symmetry);
    RUN_TEST(test_lowpass_fir_passband);
    RUN_TEST(test_lowpass_fir_stopband);
    RUN_TEST(test_lowpass_fir_near_nyquist);
    RUN_TEST(test_lowpass_fir_very_low_cutoff);
}
