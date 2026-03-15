/**
 * @file test_dtmf.c
 * @brief Tests for DTMF detection and generation.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * DTMF detection and generation
 * -----------------------------------------------------------------------*/

/** MD_dtmf_signal_length returns the correct count. */
static int test_dtmf_signal_length(void)
{
    /* 3 digits at 8 kHz, 70 ms tone, 70 ms pause:
     * tone_samples = 70 * 8000 / 1000 = 560
     * pause_samples = 560
     * total = 3*560 + 2*560 = 2800 */
    unsigned len = MD_dtmf_signal_length(3, 8000.0, 70, 70);
    if (len != 2800) return 0;

    /* Zero digits → zero samples. */
    if (MD_dtmf_signal_length(0, 8000.0, 70, 70) != 0) return 0;

    /* Single digit → just the tone, no pause. */
    unsigned one = MD_dtmf_signal_length(1, 8000.0, 70, 70);
    if (one != 560) return 0;

    return 1;
}

/** Generated DTMF tones have spectral peaks at the correct frequencies. */
static int test_dtmf_generate_frequencies(void)
{
    /* Generate digit '5' (row=770 Hz, col=1336 Hz) at sample_rate=8000. */
    const double sr = 8000.0;
    const unsigned tone_ms = 70;
    const unsigned pause_ms = 70;
    unsigned N = MD_dtmf_signal_length(1, sr, tone_ms, pause_ms);
    double *sig = malloc(N * sizeof(double));
    MD_dtmf_generate(sig, "5", sr, tone_ms, pause_ms);

    /* Compute magnitude spectrum. */
    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    /* Normalise. */
    for (unsigned k = 0; k < num_bins; k++) {
        mag[k] /= (double)N;
        if (k > 0 && k < N / 2) mag[k] *= 2.0;
    }

    /* Expected bins: 770 Hz → bin = round(770*560/8000) = 54
     *               1336 Hz → bin = round(1336*560/8000) = 94 (approx)
     * Find the two largest magnitude bins. */
    unsigned peak1 = 1, peak2 = 2;
    if (mag[peak2] > mag[peak1]) { unsigned t = peak1; peak1 = peak2; peak2 = t; }
    for (unsigned k = 3; k < num_bins; k++) {
        if (mag[k] > mag[peak1]) {
            peak2 = peak1;
            peak1 = k;
        } else if (mag[k] > mag[peak2]) {
            peak2 = k;
        }
    }

    /* Convert bins to Hz. */
    double f1 = (double)peak1 * sr / (double)N;
    double f2 = (double)peak2 * sr / (double)N;
    if (f1 > f2) { double t = f1; f1 = f2; f2 = t; }

    /* The lower peak should be near 770 Hz, upper near 1336 Hz.
     * Tolerance: ±30 Hz (one FFT bin at this resolution). */
    int ok = approx_equal(f1, 770.0, 30.0) && approx_equal(f2, 1336.0, 30.0);

    free(mag);
    free(sig);
    return ok;
}

/** Generate-then-detect round-trip recovers all digits. */
static int test_dtmf_detect_roundtrip(void)
{
    const char   *digits = "1234567890*#ABCD";
    const double  sr     = 8000.0;
    unsigned num_digits   = (unsigned)strlen(digits);
    unsigned signal_len   = MD_dtmf_signal_length(num_digits, sr, 70, 70);

    double *sig = malloc(signal_len * sizeof(double));
    MD_dtmf_generate(sig, digits, sr, 70, 70);

    MD_DTMFTone tones[32];
    unsigned n = MD_dtmf_detect(sig, signal_len, sr, tones, 32);

    int ok = (n == num_digits);
    if (ok) {
        for (unsigned i = 0; i < num_digits; i++) {
            if (tones[i].digit != digits[i]) { ok = 0; break; }
        }
    }

    free(sig);
    return ok;
}

/** Detection works with light additive noise. */
static int test_dtmf_detect_noisy(void)
{
    const char   *digits = "8675309";
    const double  sr     = 8000.0;
    unsigned num_digits   = (unsigned)strlen(digits);
    unsigned signal_len   = MD_dtmf_signal_length(num_digits, sr, 70, 70);

    double *sig   = malloc(signal_len * sizeof(double));
    double *noise = malloc(signal_len * sizeof(double));
    MD_dtmf_generate(sig, digits, sr, 70, 70);
    MD_white_noise(noise, signal_len, 0.05, 42);
    for (unsigned i = 0; i < signal_len; i++)
        sig[i] += noise[i];

    MD_DTMFTone tones[16];
    unsigned n = MD_dtmf_detect(sig, signal_len, sr, tones, 16);

    int ok = (n == num_digits);
    if (ok) {
        for (unsigned i = 0; i < num_digits; i++) {
            if (tones[i].digit != digits[i]) { ok = 0; break; }
        }
    }

    free(noise);
    free(sig);
    return ok;
}

/** Round-trip works at 16 kHz sample rate. */
static int test_dtmf_detect_16khz(void)
{
    const char   *digits = "159*0#";
    const double  sr     = 16000.0;
    unsigned num_digits   = (unsigned)strlen(digits);
    unsigned signal_len   = MD_dtmf_signal_length(num_digits, sr, 80, 80);

    double *sig = malloc(signal_len * sizeof(double));
    MD_dtmf_generate(sig, digits, sr, 80, 80);

    MD_DTMFTone tones[16];
    unsigned n = MD_dtmf_detect(sig, signal_len, sr, tones, 16);

    int ok = (n == num_digits);
    if (ok) {
        for (unsigned i = 0; i < num_digits; i++) {
            if (tones[i].digit != digits[i]) { ok = 0; break; }
        }
    }

    free(sig);
    return ok;
}

/** Tone timestamps are monotonically increasing. */
static int test_dtmf_detect_timestamps(void)
{
    const char   *digits = "159#";
    const double  sr     = 8000.0;
    unsigned signal_len   = MD_dtmf_signal_length(4, sr, 80, 60);

    double *sig = malloc(signal_len * sizeof(double));
    MD_dtmf_generate(sig, digits, sr, 80, 60);

    MD_DTMFTone tones[8];
    unsigned n = MD_dtmf_detect(sig, signal_len, sr, tones, 8);

    int ok = (n == 4);
    for (unsigned i = 0; ok && i < n; i++) {
        ok &= (tones[i].end_s > tones[i].start_s);
        if (i > 0)
            ok &= (tones[i].start_s > tones[i - 1].end_s);
    }

    free(sig);
    return ok;
}

static int test_dtmf_detect_q24_minimum(void)
{
    /* Generate at the exact Q.24 minimums: 40 ms tone, 40 ms pause. */
    const char   *digits = "19A#";
    const double  sr     = 8000.0;
    unsigned signal_len  = MD_dtmf_signal_length(4, sr, 40, 40);

    double *sig = malloc(signal_len * sizeof(double));
    MD_dtmf_generate(sig, digits, sr, 40, 40);

    MD_DTMFTone tones[8];
    unsigned n = MD_dtmf_detect(sig, signal_len, sr, tones, 8);

    int ok = (n == 4);
    for (unsigned i = 0; ok && i < n; i++)
        ok &= (tones[i].digit == digits[i]);

    free(sig);
    return ok;
}

/** Silent signal should produce zero detected digits. */
static int test_dtmf_detect_silence(void)
{
    unsigned N = 8000;  /* 1 second at 8 kHz */
    double *sig = calloc(N, sizeof(double));
    MD_DTMFTone tones[16];
    unsigned n = MD_dtmf_detect(sig, N, 8000.0, tones, 16);
    free(sig);
    return (n == 0);
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_dtmf_tests(void)
{
    printf("\n--- DTMF detection and generation ---\n");
    RUN_TEST(test_dtmf_signal_length);
    RUN_TEST(test_dtmf_generate_frequencies);
    RUN_TEST(test_dtmf_detect_roundtrip);
    RUN_TEST(test_dtmf_detect_noisy);
    RUN_TEST(test_dtmf_detect_16khz);
    RUN_TEST(test_dtmf_detect_timestamps);
    RUN_TEST(test_dtmf_detect_q24_minimum);
    RUN_TEST(test_dtmf_detect_silence);
}
