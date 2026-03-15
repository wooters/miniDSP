/**
 * @file test_resample.c
 * @brief Tests for MD_resample and MD_resample_output_len.
 */

#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * MD_resample_output_len tests
 * -----------------------------------------------------------------------*/

static int test_resample_output_len_upsample(void)
{
    return MD_resample_output_len(44100, 44100.0, 48000.0) == 48000;
}

static int test_resample_output_len_downsample(void)
{
    return MD_resample_output_len(48000, 48000.0, 44100.0) == 44100;
}

static int test_resample_output_len_noninteger(void)
{
    unsigned expected = (unsigned)ceil(1000.0 * 48000.0 / 44100.0);
    return MD_resample_output_len(1000, 44100.0, 48000.0) == expected;
}

/** Same rate should return same length. */
static int test_resample_output_len_same_rate(void)
{
    return MD_resample_output_len(1000, 44100.0, 44100.0) == 1000;
}

/* -----------------------------------------------------------------------
 * MD_resample tests
 * -----------------------------------------------------------------------*/

static int test_resample_identity(void)
{
    unsigned N = 1024;
    double *in = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));
    MD_sine_wave(in, N, 1.0, 100.0, 8000.0);

    unsigned n = MD_resample(in, N, out, N,
                             8000.0, 8000.0, 32, 10.0);
    if (n != N) { free(in); free(out); return 0; }

    for (unsigned i = 0; i < N; i++) {
        if (!approx_equal(in[i], out[i], 1e-6)) {
            free(in); free(out); return 0;
        }
    }
    free(in); free(out);
    return 1;
}

static int test_resample_dc_preservation(void)
{
    unsigned N_in = 1000;
    unsigned N_out = MD_resample_output_len(N_in, 44100.0, 48000.0);
    double *in = malloc(N_in * sizeof(double));
    double *out = malloc(N_out * sizeof(double));

    for (unsigned i = 0; i < N_in; i++) in[i] = 1.0;

    unsigned n = MD_resample(in, N_in, out, N_out,
                             44100.0, 48000.0, 32, 10.0);
    /* Interior samples (away from boundary effects) should be ~1.0.
     * The filter has 2*num_zero_crossings = 64 taps, so boundary effects
     * extend about that many samples from each edge. Use generous margin. */
    int ok = 1;
    unsigned margin = 128;
    for (unsigned i = margin; i + margin < n; i++) {
        if (!approx_equal(out[i], 1.0, 1e-4)) { ok = 0; break; }
    }
    free(in); free(out);
    return ok;
}

static int test_resample_sine_frequency(void)
{
    /* 440 Hz sine at 44100 Hz → resample to 48000 Hz → measure F0 */
    unsigned N_in = 44100;  /* 1 second */
    unsigned N_out = MD_resample_output_len(N_in, 44100.0, 48000.0);
    double *in = malloc(N_in * sizeof(double));
    double *out = malloc(N_out * sizeof(double));

    MD_sine_wave(in, N_in, 1.0, 440.0, 44100.0);

    unsigned n = MD_resample(in, N_in, out, N_out,
                             44100.0, 48000.0, 32, 10.0);

    /* Measure F0 of the output */
    double f0 = MD_f0_autocorrelation(out, n, 48000.0, 100.0, 2000.0);

    free(in); free(out);
    return fabs(f0 - 440.0) < 1.0;
}

static int test_resample_output_length_matches(void)
{
    /* Test multiple rate pairs */
    struct { double in_rate; double out_rate; } pairs[] = {
        {44100.0, 48000.0},
        {48000.0, 44100.0},
        {16000.0,  8000.0},
        { 8000.0, 16000.0},
        {44100.0, 22050.0},
    };
    unsigned N_in = 1000;

    for (unsigned p = 0; p < 5; p++) {
        unsigned expected = MD_resample_output_len(N_in, pairs[p].in_rate,
                                                   pairs[p].out_rate);
        double *in = malloc(N_in * sizeof(double));
        double *out = malloc(expected * sizeof(double));
        MD_sine_wave(in, N_in, 1.0, 100.0, pairs[p].in_rate);

        unsigned actual = MD_resample(in, N_in, out, expected,
                                      pairs[p].in_rate, pairs[p].out_rate,
                                      16, 10.0);
        free(in); free(out);
        if (actual != expected) return 0;
    }
    return 1;
}

static int test_resample_energy_preservation(void)
{
    /* White noise resampled should have similar RMS */
    unsigned N_in = 8000;
    unsigned N_out = MD_resample_output_len(N_in, 44100.0, 48000.0);
    double *in = malloc(N_in * sizeof(double));
    double *out = malloc(N_out * sizeof(double));

    MD_white_noise(in, N_in, 1.0, 42);

    unsigned n = MD_resample(in, N_in, out, N_out,
                             44100.0, 48000.0, 32, 10.0);

    double rms_in = MD_rms(in, N_in);
    /* Skip boundary samples for RMS measurement */
    unsigned margin = 128;
    double rms_out = MD_rms(out + margin, n - 2 * margin);
    double db_diff = 20.0 * log10(rms_out / rms_in);

    free(in); free(out);
    return fabs(db_diff) < 0.5;
}

static int test_resample_antialiasing(void)
{
    /* 20 kHz sine at 48 kHz, downsample to 16 kHz.
     * The 20 kHz component is above the 8 kHz Nyquist of the output,
     * so it should be suppressed. */
    unsigned N_in = 4800;  /* 0.1 seconds at 48k */
    unsigned N_out = MD_resample_output_len(N_in, 48000.0, 16000.0);
    double *in = malloc(N_in * sizeof(double));
    double *out = malloc(N_out * sizeof(double));

    MD_sine_wave(in, N_in, 1.0, 20000.0, 48000.0);

    unsigned n = MD_resample(in, N_in, out, N_out,
                             48000.0, 16000.0, 32, 10.0);

    /* Output should be near-silent (aliased frequency suppressed) */
    unsigned margin = 64;
    double rms_out = MD_rms(out + margin, n - 2 * margin);
    double attenuation_db = 20.0 * log10(fmax(rms_out, 1e-15));

    free(in); free(out);
    return attenuation_db < -60.0;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_resample_tests(void)
{
    printf("\n--- MD_resample_output_len ---\n");
    RUN_TEST(test_resample_output_len_upsample);
    RUN_TEST(test_resample_output_len_downsample);
    RUN_TEST(test_resample_output_len_noninteger);
    RUN_TEST(test_resample_output_len_same_rate);

    printf("\n--- MD_resample ---\n");
    RUN_TEST(test_resample_identity);
    RUN_TEST(test_resample_dc_preservation);
    RUN_TEST(test_resample_sine_frequency);
    RUN_TEST(test_resample_output_length_matches);
    RUN_TEST(test_resample_energy_preservation);
    RUN_TEST(test_resample_antialiasing);
}
