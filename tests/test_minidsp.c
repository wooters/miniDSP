/**
 * @file test_minidsp.c
 * @brief Comprehensive test suite for the miniDSP library.
 *
 * This file tests every public function in minidsp.h and biquad.h
 * using known mathematical properties to verify correctness.
 *
 * Each test function returns 1 on success and 0 on failure.
 * The main() function runs all tests and prints a summary.
 *
 * How to compile (from the tests/ directory):
 *   make test_minidsp
 *
 * How to run:
 *   ./test_minidsp
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <unistd.h>
#include "minidsp.h"
#include "biquad.h"
#include "fileio.h"

/* -----------------------------------------------------------------------
 * Test infrastructure
 *
 * These macros and counters keep track of how many tests pass and fail.
 * -----------------------------------------------------------------------*/

static int tests_run    = 0;
static int tests_passed = 0;
static int tests_failed = 0;

/** Compare two doubles with a tolerance.  Returns 1 if close enough. */
static int approx_equal(double a, double b, double tolerance)
{
    return fabs(a - b) <= tolerance;
}

/** Run a single test and print the result. */
#define RUN_TEST(test_func)                                       \
    do {                                                          \
        tests_run++;                                              \
        printf("  %-50s", #test_func);                            \
        if (test_func()) {                                        \
            tests_passed++;                                       \
            printf("[PASS]\n");                                   \
        } else {                                                  \
            tests_failed++;                                       \
            printf("[FAIL]\n");                                   \
        }                                                         \
    } while (0)

/* -----------------------------------------------------------------------
 * Helper: create a delayed copy of a signal using circular rotation.
 *
 * This mimics what happens when a sound reaches one microphone later
 * than another.  A positive 'delay' shifts the signal to the right
 * (later in time).
 * -----------------------------------------------------------------------*/
static void delay_signal(const double *in, double *out,
                         unsigned n, int delay)
{
    for (unsigned i = 0; i < n; i++) {
        /* Use modular arithmetic with unsigned wraparound.
         * For negative delays, the unsigned conversion + modulus
         * produces the correct circular shift. */
        unsigned j = (i + (unsigned)delay) % n;
        out[j] = in[i];
    }
}

/* -----------------------------------------------------------------------
 * Tests for MD_dot()
 * -----------------------------------------------------------------------*/

/** Dot product of orthogonal vectors should be zero. */
static int test_dot_orthogonal(void)
{
    double a[] = {1.0, 0.0, 0.0};
    double b[] = {0.0, 1.0, 0.0};
    return approx_equal(MD_dot(a, b, 3), 0.0, 1e-15);
}

/** Dot product of a vector with itself equals the sum of squares. */
static int test_dot_self(void)
{
    double a[] = {1.0, 2.0, 3.0};
    /* 1^2 + 2^2 + 3^2 = 14 */
    return approx_equal(MD_dot(a, a, 3), 14.0, 1e-15);
}

/** Dot product of known vectors. */
static int test_dot_known(void)
{
    double a[] = {1.0, 2.0, 3.0, 4.0};
    double b[] = {5.0, 6.0, 7.0, 8.0};
    /* 1*5 + 2*6 + 3*7 + 4*8 = 5+12+21+32 = 70 */
    return approx_equal(MD_dot(a, b, 4), 70.0, 1e-15);
}

/** Dot product with length 1 just returns a*b. */
static int test_dot_single(void)
{
    double a[] = {3.5};
    double b[] = {2.0};
    return approx_equal(MD_dot(a, b, 1), 7.0, 1e-15);
}

/** Dot product with all zeros should be zero. */
static int test_dot_zeros(void)
{
    double a[] = {0.0, 0.0, 0.0};
    double b[] = {1.0, 2.0, 3.0};
    return approx_equal(MD_dot(a, b, 3), 0.0, 1e-15);
}

/* -----------------------------------------------------------------------
 * Tests for MD_energy()
 * -----------------------------------------------------------------------*/

/** Energy of a zero signal is zero. */
static int test_energy_zero(void)
{
    double a[] = {0.0, 0.0, 0.0};
    return approx_equal(MD_energy(a, 3), 0.0, 1e-15);
}

/** Energy of a constant signal: N * c^2. */
static int test_energy_constant(void)
{
    double a[] = {3.0, 3.0, 3.0, 3.0};
    /* 4 * 3^2 = 36 */
    return approx_equal(MD_energy(a, 4), 36.0, 1e-15);
}

/** Energy of a known signal. */
static int test_energy_known(void)
{
    double a[] = {1.0, -2.0, 3.0};
    /* 1 + 4 + 9 = 14 */
    return approx_equal(MD_energy(a, 3), 14.0, 1e-15);
}

/** Energy of a single sample. */
static int test_energy_single(void)
{
    double a[] = {5.0};
    return approx_equal(MD_energy(a, 1), 25.0, 1e-15);
}

/* -----------------------------------------------------------------------
 * Tests for MD_power()
 * -----------------------------------------------------------------------*/

/** Power = energy / N. */
static int test_power_known(void)
{
    double a[] = {1.0, -2.0, 3.0};
    /* energy = 14, N = 3, power = 14/3 */
    return approx_equal(MD_power(a, 3), 14.0 / 3.0, 1e-15);
}

/** Power of a unit-amplitude sine wave over full periods approaches 0.5. */
static int test_power_sine(void)
{
    /* Generate exactly one full cycle of sin(x) with many samples */
    unsigned N = 10000;
    double *sig = malloc(N * sizeof(double));
    for (unsigned i = 0; i < N; i++) {
        sig[i] = sin(2.0 * M_PI * (double)i / (double)N);
    }
    /* Power of sin(x) over one full cycle = 0.5 */
    double p = MD_power(sig, N);
    free(sig);
    return approx_equal(p, 0.5, 1e-4);
}

/* -----------------------------------------------------------------------
 * Tests for MD_power_db()
 * -----------------------------------------------------------------------*/

/** Power in dB for a known signal. */
static int test_power_db_known(void)
{
    /* Create a signal with power = 1.0 -> 0 dB */
    double a[] = {1.0, -1.0, 1.0, -1.0};
    /* energy = 4, N = 4, power = 1.0, dB = 0.0 */
    return approx_equal(MD_power_db(a, 4), 0.0, 1e-10);
}

/** A very quiet signal should be at the dB floor (power gets clamped to 1e-10). */
static int test_power_db_quiet(void)
{
    double a[] = {1e-6, 1e-6};
    /* power = 1e-12, but the floor of 1e-10 kicks in.
     * dB = 10*log10(1e-10) = -100.0 */
    double db = MD_power_db(a, 2);
    return (db <= -100.0);
}

/** Silent signal should hit the floor (not -inf). */
static int test_power_db_floor(void)
{
    double a[] = {0.0, 0.0, 0.0};
    double db = MD_power_db(a, 3);
    /* The floor is 1e-10, so dB = 10*log10(1e-10) = -100 */
    return approx_equal(db, -100.0, 1e-10);
}

/* -----------------------------------------------------------------------
 * Tests for MD_scale() and MD_scale_vec()
 * -----------------------------------------------------------------------*/

/** Scale a midpoint value. */
static int test_scale_midpoint(void)
{
    /* 5 is the midpoint of [0,10], should map to midpoint of [0,100] = 50 */
    return approx_equal(MD_scale(5.0, 0.0, 10.0, 0.0, 100.0), 50.0, 1e-15);
}

/** Scale endpoints. */
static int test_scale_endpoints(void)
{
    int ok = 1;
    ok &= approx_equal(MD_scale(0.0, 0.0, 10.0, 0.0, 100.0), 0.0, 1e-15);
    ok &= approx_equal(MD_scale(10.0, 0.0, 10.0, 0.0, 100.0), 100.0, 1e-15);
    return ok;
}

/** Scale a vector. */
static int test_scale_vec(void)
{
    double in[]  = {0.0, 5.0, 10.0};
    double out[3];
    MD_scale_vec(in, out, 3, 0.0, 10.0, -1.0, 1.0);
    int ok = 1;
    ok &= approx_equal(out[0], -1.0, 1e-15);
    ok &= approx_equal(out[1],  0.0, 1e-15);
    ok &= approx_equal(out[2],  1.0, 1e-15);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_fit_within_range()
 * -----------------------------------------------------------------------*/

/** Values already in range should be copied unchanged. */
static int test_fit_in_range_no_change(void)
{
    double in[]  = {0.1, 0.5, 0.9};
    double out[3];
    MD_fit_within_range(in, out, 3, 0.0, 1.0);
    int ok = 1;
    ok &= approx_equal(out[0], 0.1, 1e-15);
    ok &= approx_equal(out[1], 0.5, 1e-15);
    ok &= approx_equal(out[2], 0.9, 1e-15);
    return ok;
}

/** Values outside range should be rescaled. */
static int test_fit_in_range_rescale(void)
{
    double in[]  = {-10.0, 0.0, 10.0};
    double out[3];
    MD_fit_within_range(in, out, 3, -1.0, 1.0);
    int ok = 1;
    ok &= approx_equal(out[0], -1.0, 1e-15);
    ok &= approx_equal(out[1],  0.0, 1e-15);
    ok &= approx_equal(out[2],  1.0, 1e-15);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_adjust_dblevel()
 * -----------------------------------------------------------------------*/

/** After adjustment, the output power should match the target dB. */
static int test_adjust_dblevel(void)
{
    unsigned N = 1000;
    double *in  = malloc(N * sizeof(double));
    double *out = malloc(N * sizeof(double));

    /* Create a low-level signal */
    for (unsigned i = 0; i < N; i++) {
        in[i] = 0.01 * sin(2.0 * M_PI * (double)i / 100.0);
    }

    double target_db = -10.0;
    MD_adjust_dblevel(in, out, N, target_db);

    double actual_db = MD_power_db(out, N);

    free(out);
    free(in);

    /* The output dB might not exactly match if clipping occurred,
     * but for a small signal boosted to -10 dB it should be close. */
    return approx_equal(actual_db, target_db, 0.5);
}

/* -----------------------------------------------------------------------
 * Tests for MD_entropy()
 * -----------------------------------------------------------------------*/

/** Uniform distribution should have maximum entropy (close to 1.0). */
static int test_entropy_uniform(void)
{
    double a[] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    double e = MD_entropy(a, 8, true);
    return approx_equal(e, 1.0, 1e-10);
}

/** Spike distribution (all energy in one bin) should have low entropy. */
static int test_entropy_spike(void)
{
    double a[] = {0.0, 0.0, 1000.0, 0.0, 0.0};
    double e = MD_entropy(a, 5, true);
    return approx_equal(e, 0.0, 1e-10);
}

/** All-zeros should return 0.0 (no meaningful distribution). */
static int test_entropy_zeros(void)
{
    double a[] = {0.0, 0.0, 0.0};
    double e = MD_entropy(a, 3, true);
    return approx_equal(e, 0.0, 1e-15);
}

/** Single element should return 0.0. */
static int test_entropy_single(void)
{
    double a[] = {42.0};
    double e = MD_entropy(a, 1, true);
    return approx_equal(e, 0.0, 1e-15);
}

/** Non-clip mode squares values (handles negatives). */
static int test_entropy_no_clip(void)
{
    /* Symmetric distribution: squared values should give uniform-ish distribution */
    double a[] = {-1.0, 1.0, -1.0, 1.0};
    double e = MD_entropy(a, 4, false);
    /* All squared values are 1.0, so this is a uniform distribution -> entropy = 1.0 */
    return approx_equal(e, 1.0, 1e-10);
}

/** Entropy with clip mode should ignore negatives. */
static int test_entropy_clip(void)
{
    double a[] = {-5.0, -3.0, 0.0, 1.0};
    double e = MD_entropy(a, 4, true);
    /* Only a[3]=1.0 contributes.  p=1.0, -p*log2(p) = 0.  Normalized: 0. */
    return approx_equal(e, 0.0, 1e-10);
}

/* -----------------------------------------------------------------------
 * Tests for MD_rms()
 * -----------------------------------------------------------------------*/

/** RMS of a unit-amplitude sine wave should be 1/sqrt(2). */
static int test_rms_sine(void)
{
    unsigned N = 8000;
    double sig[8000];
    /* Use integer number of complete cycles to avoid edge effects */
    MD_sine_wave(sig, N, 1.0, 100.0, 8000.0);
    double rms = MD_rms(sig, N);
    return approx_equal(rms, 1.0 / sqrt(2.0), 1e-6);
}

/** RMS of a DC signal of value c should be |c|. */
static int test_rms_dc(void)
{
    double sig[100];
    for (unsigned i = 0; i < 100; i++) sig[i] = 3.0;
    return approx_equal(MD_rms(sig, 100), 3.0, 1e-15);
}

/** RMS of silence should be 0. */
static int test_rms_silence(void)
{
    double sig[64];
    memset(sig, 0, sizeof(sig));
    return approx_equal(MD_rms(sig, 64), 0.0, 1e-15);
}

/** RMS should equal sqrt(MD_power()). */
static int test_rms_matches_sqrt_power(void)
{
    double sig[] = {1.0, -2.0, 3.0, -4.0, 5.0};
    double rms = MD_rms(sig, 5);
    double expected = sqrt(MD_power(sig, 5));
    return approx_equal(rms, expected, 1e-15);
}

/* -----------------------------------------------------------------------
 * Tests for MD_zero_crossing_rate()
 * -----------------------------------------------------------------------*/

/** ZCR of a sine wave should be approximately 2*freq/sample_rate. */
static int test_zcr_sine(void)
{
    unsigned N = 8000;
    double sig[8000];
    double freq = 500.0, sr = 8000.0;
    MD_sine_wave(sig, N, 1.0, freq, sr);
    double zcr = MD_zero_crossing_rate(sig, N);
    double expected = 2.0 * freq / sr;  /* 0.125 */
    return approx_equal(zcr, expected, 0.005);
}

/** ZCR of a constant signal should be 0. */
static int test_zcr_constant(void)
{
    double sig[100];
    for (unsigned i = 0; i < 100; i++) sig[i] = 5.0;
    return approx_equal(MD_zero_crossing_rate(sig, 100), 0.0, 1e-15);
}

/** ZCR of a perfectly alternating signal should be 1.0. */
static int test_zcr_alternating(void)
{
    double sig[100];
    for (unsigned i = 0; i < 100; i++) sig[i] = (i % 2 == 0) ? 1.0 : -1.0;
    return approx_equal(MD_zero_crossing_rate(sig, 100), 1.0, 1e-15);
}

/** Noise should have higher ZCR than a low-frequency sine. */
static int test_zcr_noise_higher_than_sine(void)
{
    unsigned N = 4096;
    double sine[4096], noise[4096];
    MD_sine_wave(sine, N, 1.0, 100.0, 8000.0);
    MD_white_noise(noise, N, 1.0, 42);
    return MD_zero_crossing_rate(noise, N) > MD_zero_crossing_rate(sine, N);
}

/** Zeros in the signal should be treated as non-negative. */
static int test_zcr_zeros_handling(void)
{
    /* {1, 0, -1, 0, 1}: sign changes at indices 1->2 and 2->3 */
    double sig[] = {1.0, 0.0, -1.0, 0.0, 1.0};
    double zcr = MD_zero_crossing_rate(sig, 5);
    /* 0 is non-negative, so transitions: 0→-1 (cross), -1→0 (cross) = 2/4 = 0.5 */
    return approx_equal(zcr, 0.5, 1e-15);
}

/* -----------------------------------------------------------------------
 * Tests for MD_autocorrelation()
 * -----------------------------------------------------------------------*/

/** Autocorrelation at lag 0 should always be 1.0. */
static int test_acf_lag0(void)
{
    double sig[] = {1.0, -2.0, 3.0, -4.0, 5.0};
    double acf[3];
    MD_autocorrelation(sig, 5, acf, 3);
    return approx_equal(acf[0], 1.0, 1e-15);
}

/** Sine autocorrelation should peak at the period. */
static int test_acf_sine_period(void)
{
    unsigned N = 2048;
    double sig[2048];
    /* 100 Hz at 1000 Hz sample rate -> period = 10 samples */
    MD_sine_wave(sig, N, 1.0, 100.0, 1000.0);
    unsigned max_lag = 50;
    double acf[50];
    MD_autocorrelation(sig, N, acf, max_lag);
    /* acf[10] should be close to 1.0 (one full period) */
    return approx_equal(acf[10], 1.0, 0.01);
}

/** Noise autocorrelation should decay quickly. */
static int test_acf_noise_decay(void)
{
    unsigned N = 4096;
    double *sig = malloc(N * sizeof(double));
    MD_white_noise(sig, N, 1.0, 42);
    unsigned max_lag = 100;
    double *acf = malloc(max_lag * sizeof(double));
    MD_autocorrelation(sig, N, acf, max_lag);
    /* acf[0] = 1.0, acf at larger lags should be much smaller */
    int ok = (fabs(acf[10]) < 0.15);
    free(acf);
    free(sig);
    return ok;
}

/** All autocorrelation values should have magnitude <= 1.0. */
static int test_acf_bounded(void)
{
    double sig[] = {1.0, -3.0, 2.0, -1.0, 4.0, -2.0, 3.0, -5.0};
    double acf[5];
    MD_autocorrelation(sig, 8, acf, 5);
    for (unsigned i = 0; i < 5; i++) {
        if (fabs(acf[i]) > 1.0 + 1e-10) return 0;
    }
    return 1;
}

/** Autocorrelation of silence should be all zeros. */
static int test_acf_silence(void)
{
    double sig[32];
    memset(sig, 0, sizeof(sig));
    double acf[10];
    MD_autocorrelation(sig, 32, acf, 10);
    for (unsigned i = 0; i < 10; i++) {
        if (fabs(acf[i]) > 1e-15) return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * Tests for MD_peak_detect()
 * -----------------------------------------------------------------------*/

/** Detect known peaks in a simple signal. */
static int test_peaks_known(void)
{
    double sig[] = {0, 1, 3, 1, 0, 2, 5, 2, 0};
    unsigned peaks[9], n;
    MD_peak_detect(sig, 9, 0.0, 1, peaks, &n);
    return n == 2 && peaks[0] == 2 && peaks[1] == 6;
}

/** Threshold should filter out small peaks. */
static int test_peaks_threshold(void)
{
    double sig[] = {0, 1, 3, 1, 0, 2, 5, 2, 0};
    unsigned peaks[9], n;
    MD_peak_detect(sig, 9, 4.0, 1, peaks, &n);
    /* Only the peak at index 6 (value 5) is above threshold 4 */
    return n == 1 && peaks[0] == 6;
}

/** min_distance should suppress nearby peaks. */
static int test_peaks_min_distance(void)
{
    double sig[] = {0, 3, 0, 2, 0, 4, 0};
    unsigned peaks[7], n;
    MD_peak_detect(sig, 7, 0.0, 3, peaks, &n);
    /* Peak at 1 (value 3) accepted, peak at 3 (value 2) is within distance 3, skipped.
       Peak at 5 (value 4) is at distance 4 from peak 1, accepted. */
    return n == 2 && peaks[0] == 1 && peaks[1] == 5;
}

/** Signal with no peaks should return zero count. */
static int test_peaks_none(void)
{
    double sig[] = {1, 2, 3, 4, 5};
    unsigned peaks[5], n;
    MD_peak_detect(sig, 5, 0.0, 1, peaks, &n);
    /* Monotonically increasing -> no local maxima */
    return n == 0;
}

/** Flat signal should produce zero peaks (equal neighbours). */
static int test_peaks_flat_signal(void)
{
    double sig[] = {3, 3, 3, 3, 3};
    unsigned peaks[5], n;
    MD_peak_detect(sig, 5, 0.0, 1, peaks, &n);
    return n == 0;
}

/** Single-sample signal should produce zero peaks. */
static int test_peaks_single_sample(void)
{
    double sig[] = {5.0};
    unsigned peaks[1], n;
    MD_peak_detect(sig, 1, 0.0, 1, peaks, &n);
    return n == 0;
}

/* -----------------------------------------------------------------------
 * Tests for pitch detection (MD_f0_autocorrelation, MD_f0_fft)
 * -----------------------------------------------------------------------*/

/** Autocorrelation F0 should recover a clean sine tone. */
static int test_f0_acf_clean_sine(void)
{
    const unsigned N = 4096;
    const double fs = 16000.0;
    const double f0 = 220.0;
    double sig[N];

    MD_sine_wave(sig, N, 1.0, f0, fs);
    double est = MD_f0_autocorrelation(sig, N, fs, 80.0, 400.0);
    return approx_equal(est, f0, 1.0);
}

/** Autocorrelation should prefer the fundamental on harmonic-rich signals. */
static int test_f0_acf_harmonic_signal(void)
{
    const unsigned N = 4096;
    const double fs = 8000.0;
    const double f0 = 200.0;
    double sig[N];

    MD_square_wave(sig, N, 1.0, f0, fs);
    double est = MD_f0_autocorrelation(sig, N, fs, 80.0, 400.0);
    return approx_equal(est, f0, 2.0);
}

/** Autocorrelation F0 should remain close on a noisy sine. */
static int test_f0_acf_noisy_sine(void)
{
    const unsigned N = 4096;
    const double fs = 16000.0;
    const double f0 = 220.0;
    double sig[N], noise[N];

    MD_sine_wave(sig, N, 0.8, f0, fs);
    MD_white_noise(noise, N, 0.2, 42);
    for (unsigned i = 0; i < N; i++) {
        sig[i] += noise[i];
    }

    double est = MD_f0_autocorrelation(sig, N, fs, 80.0, 400.0);
    return approx_equal(est, f0, 5.0);
}

/** Autocorrelation F0 should return 0 for silence. */
static int test_f0_acf_silence(void)
{
    double sig[2048] = {0};
    double est = MD_f0_autocorrelation(sig, 2048, 16000.0, 80.0, 400.0);
    return approx_equal(est, 0.0, 1e-15);
}

/** Autocorrelation F0 should return 0 when the true F0 is outside range. */
static int test_f0_acf_out_of_range(void)
{
    const unsigned N = 4096;
    double sig[N];
    MD_sine_wave(sig, N, 1.0, 220.0, 16000.0);
    double est = MD_f0_autocorrelation(sig, N, 16000.0, 300.0, 500.0);
    return approx_equal(est, 0.0, 1e-15);
}

/** Autocorrelation F0 should handle lag-bound clamping at both edges. */
static int test_f0_acf_lag_edge_mapping(void)
{
    const unsigned N = 1024;
    const double fs = 8000.0;
    const double f0 = 500.0;
    double sig[N];

    MD_sine_wave(sig, N, 1.0, f0, fs);
    /* min/max intentionally exceed practical bounds to trigger clamping:
     * lag_min -> 1, lag_max -> N-1. */
    double est = MD_f0_autocorrelation(sig, N, fs, 1.0, 20000.0);
    return approx_equal(est, f0, 2.0);
}

/** FFT F0 should recover a clean sine tone. */
static int test_f0_fft_clean_sine(void)
{
    const unsigned N = 4096;
    const double fs = 16000.0;
    const double f0 = 220.0;
    double sig[N];

    MD_sine_wave(sig, N, 1.0, f0, fs);
    double est = MD_f0_fft(sig, N, fs, 80.0, 400.0);
    return approx_equal(est, f0, 2.0);
}

/** FFT F0 should pick the dominant tone when two tones are present. */
static int test_f0_fft_two_tone_dominance(void)
{
    const unsigned N = 4096;
    const double fs = 16000.0;
    double tone_a[N], tone_b[N], mix[N];

    MD_sine_wave(tone_a, N, 1.0, 220.0, fs);
    MD_sine_wave(tone_b, N, 0.35, 330.0, fs);
    for (unsigned i = 0; i < N; i++) {
        mix[i] = tone_a[i] + tone_b[i];
    }

    double est = MD_f0_fft(mix, N, fs, 80.0, 400.0);
    return approx_equal(est, 220.0, 2.0);
}

/** FFT F0 should return 0 for silence. */
static int test_f0_fft_silence(void)
{
    double sig[2048] = {0};
    double est = MD_f0_fft(sig, 2048, 16000.0, 80.0, 400.0);
    return approx_equal(est, 0.0, 1e-15);
}

/** FFT F0 should return 0 when the true F0 is outside range. */
static int test_f0_fft_out_of_range(void)
{
    const unsigned N = 4096;
    double sig[N];
    MD_sine_wave(sig, N, 1.0, 220.0, 16000.0);
    double est = MD_f0_fft(sig, N, 16000.0, 300.0, 500.0);
    return approx_equal(est, 0.0, 1e-15);
}

/** FFT F0 should handle bin-bound clamping at DC/Nyquist limits. */
static int test_f0_fft_bin_edge_mapping(void)
{
    const unsigned N = 2048;
    const double fs = 8192.0;
    const double f0 = 84.0;
    double sig[N];

    MD_sine_wave(sig, N, 1.0, f0, fs);
    double est = MD_f0_fft(sig, N, fs, 0.5, 50000.0);
    return approx_equal(est, f0, 3.0);
}

/** FFT F0 should work across consecutive calls with different frame sizes. */
static int test_f0_fft_plan_recache(void)
{
    int ok = 1;

    {
        const unsigned N = 4096;
        double sig[N];
        MD_sine_wave(sig, N, 1.0, 220.0, 16000.0);
        double est = MD_f0_fft(sig, N, 16000.0, 80.0, 400.0);
        ok &= approx_equal(est, 220.0, 2.0);
    }

    {
        const unsigned N = 1024;
        double sig[N];
        MD_sine_wave(sig, N, 1.0, 300.0, 8000.0);
        double est = MD_f0_fft(sig, N, 8000.0, 80.0, 400.0);
        ok &= approx_equal(est, 300.0, 4.0);
    }

    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_mix()
 * -----------------------------------------------------------------------*/

/** Equal weights of 0.5 should produce the average. */
static int test_mix_equal_weights(void)
{
    double a[] = {2.0, 4.0, 6.0};
    double b[] = {10.0, 20.0, 30.0};
    double out[3];
    MD_mix(a, b, out, 3, 0.5, 0.5);
    int ok = 1;
    ok &= approx_equal(out[0], 6.0, 1e-15);
    ok &= approx_equal(out[1], 12.0, 1e-15);
    ok &= approx_equal(out[2], 18.0, 1e-15);
    return ok;
}

/** Weight 1.0/0.0 should pass through the first signal. */
static int test_mix_passthrough(void)
{
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {10.0, 20.0, 30.0};
    double out[3];
    MD_mix(a, b, out, 3, 1.0, 0.0);
    int ok = 1;
    ok &= approx_equal(out[0], 1.0, 1e-15);
    ok &= approx_equal(out[1], 2.0, 1e-15);
    ok &= approx_equal(out[2], 3.0, 1e-15);
    return ok;
}

/** Mixing should satisfy the expected energy relationship for orthogonal signals. */
static int test_mix_energy(void)
{
    unsigned N = 8000;
    double sine[8000], noise[8000], mix[8000];
    /* Sine at an exact integer bin and noise are approximately orthogonal */
    MD_sine_wave(sine, N, 1.0, 100.0, 8000.0);
    MD_white_noise(noise, N, 0.5, 123);
    MD_mix(sine, noise, mix, N, 1.0, 1.0);
    /* Energy of mix should be approximately E(sine) + E(noise) for uncorrelated signals */
    double e_sine = MD_energy(sine, N);
    double e_noise = MD_energy(noise, N);
    double e_mix = MD_energy(mix, N);
    return approx_equal(e_mix, e_sine + e_noise, e_mix * 0.1);
}

/** In-place mixing: output aliases input a. */
static int test_mix_inplace(void)
{
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {4.0, 5.0, 6.0};
    MD_mix(a, b, a, 3, 0.5, 0.5);
    int ok = 1;
    ok &= approx_equal(a[0], 2.5, 1e-15);
    ok &= approx_equal(a[1], 3.5, 1e-15);
    ok &= approx_equal(a[2], 4.5, 1e-15);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for simple effects
 * -----------------------------------------------------------------------*/

/** Delay echo with wet=0 should pass through dry input exactly. */
static int test_delay_echo_dry_passthrough(void)
{
    double in[] = {0.1, -0.2, 0.3, -0.4, 0.5};
    double out[5];
    MD_delay_echo(in, out, 5, 2, 0.45, 1.0, 0.0);

    int ok = 1;
    for (unsigned i = 0; i < 5; i++) {
        ok &= approx_equal(out[i], in[i], 1e-15);
    }
    return ok;
}

/** Delay echo impulse response should repeat every delay with geometric decay. */
static int test_delay_echo_impulse_decay(void)
{
    const unsigned N = 20;
    double in[N];
    double out[N];
    memset(in, 0, sizeof(in));
    in[0] = 1.0;

    MD_delay_echo(in, out, N, 4, 0.5, 0.0, 1.0);

    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        double expected = 0.0;
        if (i >= 4 && i % 4 == 0) {
            unsigned m = i / 4 - 1;
            expected = pow(0.5, (double)m);
        }
        ok &= approx_equal(out[i], expected, 1e-12);
    }
    return ok;
}

/** Delay echo should be in-place safe. */
static int test_delay_echo_inplace(void)
{
    double inout[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ref_out[6];
    double ref_in[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    MD_delay_echo(ref_in, ref_out, 6, 2, 0.5, 0.0, 1.0);
    MD_delay_echo(inout, inout, 6, 2, 0.5, 0.0, 1.0);

    int ok = 1;
    for (unsigned i = 0; i < 6; i++) {
        ok &= approx_equal(inout[i], ref_out[i], 1e-12);
    }
    return ok;
}

/** Tremolo depth 0 should be exact passthrough. */
static int test_tremolo_depth_zero_passthrough(void)
{
    double in[] = {0.1, -0.2, 0.3, -0.4, 0.5};
    double out[5];
    MD_tremolo(in, out, 5, 5.0, 0.0, 8000.0);

    int ok = 1;
    for (unsigned i = 0; i < 5; i++) {
        ok &= approx_equal(out[i], in[i], 1e-15);
    }
    return ok;
}

/** Tremolo gain should stay in [1-depth, 1] for positive input. */
static int test_tremolo_gain_bounds(void)
{
    const unsigned N = 4000;
    double in[N], out[N];
    for (unsigned i = 0; i < N; i++) in[i] = 1.0;

    double depth = 0.8;
    MD_tremolo(in, out, N, 5.0, depth, 8000.0);

    double lo = 1.0 - depth;
    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        ok &= (out[i] >= lo - 1e-12);
        ok &= (out[i] <= 1.0 + 1e-12);
    }
    return ok;
}

/** Comb reverb with wet=0 should pass through dry input exactly. */
static int test_comb_reverb_dry_passthrough(void)
{
    double in[] = {0.2, -0.1, 0.0, 0.4, -0.3};
    double out[5];
    MD_comb_reverb(in, out, 5, 3, 0.75, 1.0, 0.0);

    int ok = 1;
    for (unsigned i = 0; i < 5; i++) {
        ok &= approx_equal(out[i], in[i], 1e-15);
    }
    return ok;
}

/** Comb reverb impulse should produce geometric repeats at delay multiples. */
static int test_comb_reverb_impulse_decay(void)
{
    const unsigned N = 20;
    double in[N];
    double out[N];
    memset(in, 0, sizeof(in));
    in[0] = 1.0;

    MD_comb_reverb(in, out, N, 4, 0.5, 0.0, 1.0);

    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        double expected = 0.0;
        if (i % 4 == 0) {
            unsigned m = i / 4;
            expected = pow(0.5, (double)m);
        }
        ok &= approx_equal(out[i], expected, 1e-12);
    }
    return ok;
}

/** Comb reverb should be in-place safe. */
static int test_comb_reverb_inplace(void)
{
    double inout[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    double ref_out[6];
    double ref_in[] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    MD_comb_reverb(ref_in, ref_out, 6, 2, 0.5, 0.0, 1.0);
    MD_comb_reverb(inout, inout, 6, 2, 0.5, 0.0, 1.0);

    int ok = 1;
    for (unsigned i = 0; i < 6; i++) {
        ok &= approx_equal(inout[i], ref_out[i], 1e-12);
    }
    return ok;
}

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

/* -----------------------------------------------------------------------
 * Tests for window generation
 * -----------------------------------------------------------------------*/

/** Hanning window endpoints should be zero. */
static int test_hann_endpoints(void)
{
    double win[64];
    MD_Gen_Hann_Win(win, 64);
    int ok = 1;
    ok &= approx_equal(win[0], 0.0, 1e-15);
    ok &= approx_equal(win[63], 0.0, 1e-15);
    return ok;
}

/** Hanning window peak should be 1.0 at the centre. */
static int test_hann_peak(void)
{
    unsigned n = 65;  /* odd length so there's an exact center */
    double win[65];
    MD_Gen_Hann_Win(win, n);
    return approx_equal(win[32], 1.0, 1e-10);
}

/** Hanning window should be symmetric. */
static int test_hann_symmetry(void)
{
    unsigned n = 128;
    double win[128];
    MD_Gen_Hann_Win(win, n);
    int ok = 1;
    for (unsigned i = 0; i < n / 2; i++) {
        ok &= approx_equal(win[i], win[n - 1 - i], 1e-14);
    }
    return ok;
}

/** All Hanning window values should be in [0, 1]. */
static int test_hann_range(void)
{
    unsigned n = 256;
    double win[256];
    MD_Gen_Hann_Win(win, n);
    for (unsigned i = 0; i < n; i++) {
        if (win[i] < -1e-15 || win[i] > 1.0 + 1e-15) return 0;
    }
    return 1;
}

/** Hamming window endpoints should be 0.08. */
static int test_hamming_endpoints(void)
{
    double win[64];
    MD_Gen_Hamming_Win(win, 64);
    int ok = 1;
    ok &= approx_equal(win[0], 0.08, 1e-14);
    ok &= approx_equal(win[63], 0.08, 1e-14);
    return ok;
}

/** Hamming window peak should be 1.0 at the centre. */
static int test_hamming_peak(void)
{
    unsigned n = 65;  /* odd length so there's an exact center */
    double win[65];
    MD_Gen_Hamming_Win(win, n);
    return approx_equal(win[32], 1.0, 1e-12);
}

/** Hamming window should be symmetric. */
static int test_hamming_symmetry(void)
{
    unsigned n = 128;
    double win[128];
    MD_Gen_Hamming_Win(win, n);
    int ok = 1;
    for (unsigned i = 0; i < n / 2; i++) {
        ok &= approx_equal(win[i], win[n - 1 - i], 1e-14);
    }
    return ok;
}

/** Hamming window values should be in [0.08, 1]. */
static int test_hamming_range(void)
{
    unsigned n = 256;
    double win[256];
    MD_Gen_Hamming_Win(win, n);
    for (unsigned i = 0; i < n; i++) {
        if (win[i] < 0.08 - 1e-12 || win[i] > 1.0 + 1e-12) return 0;
    }
    return 1;
}

/** Blackman window endpoints should be zero. */
static int test_blackman_endpoints(void)
{
    double win[64];
    MD_Gen_Blackman_Win(win, 64);
    int ok = 1;
    ok &= approx_equal(win[0], 0.0, 1e-14);
    ok &= approx_equal(win[63], 0.0, 1e-14);
    return ok;
}

/** Blackman window should be symmetric. */
static int test_blackman_symmetry(void)
{
    unsigned n = 128;
    double win[128];
    MD_Gen_Blackman_Win(win, n);
    int ok = 1;
    for (unsigned i = 0; i < n / 2; i++) {
        ok &= approx_equal(win[i], win[n - 1 - i], 1e-14);
    }
    return ok;
}

/** Blackman window values should be in [0, 1] (allow tiny numeric noise). */
static int test_blackman_range(void)
{
    unsigned n = 256;
    double win[256];
    MD_Gen_Blackman_Win(win, n);
    for (unsigned i = 0; i < n; i++) {
        if (win[i] < -1e-12 || win[i] > 1.0 + 1e-12) return 0;
    }
    return 1;
}

/** Rectangular window should be all ones. */
static int test_rect_all_ones(void)
{
    unsigned n = 256;
    double win[256];
    MD_Gen_Rect_Win(win, n);
    for (unsigned i = 0; i < n; i++) {
        if (!approx_equal(win[i], 1.0, 1e-15)) return 0;
    }
    return 1;
}

/** n=1 should produce a single 1.0 sample for all windows. */
static int test_window_single_sample(void)
{
    double hann = 0.0, hamming = 0.0, blackman = 0.0, rect = 0.0;
    MD_Gen_Hann_Win(&hann, 1);
    MD_Gen_Hamming_Win(&hamming, 1);
    MD_Gen_Blackman_Win(&blackman, 1);
    MD_Gen_Rect_Win(&rect, 1);
    return approx_equal(hann, 1.0, 1e-15)
        && approx_equal(hamming, 1.0, 1e-15)
        && approx_equal(blackman, 1.0, 1e-15)
        && approx_equal(rect, 1.0, 1e-15);
}

/* -----------------------------------------------------------------------
 * Shared helpers for magnitude/PSD test pairs
 *
 * MD_magnitude_spectrum and MD_power_spectral_density share the same
 * (const double *, unsigned, double *) signature.  Four tests are
 * identical except for which function is called, so we parameterize
 * them with a function pointer.
 * -----------------------------------------------------------------------*/

typedef void (*spectrum_fn_t)(const double *, unsigned, double *);

/** A single-frequency sine should produce a peak at the correct bin. */
static int _test_spectrum_single_sine(spectrum_fn_t fn)
{
    unsigned N = 1024;
    double sample_rate = 1024.0;
    double freq = 100.0;

    double *sig = malloc(N * sizeof(double));
    for (unsigned i = 0; i < N; i++) {
        sig[i] = sin(2.0 * M_PI * freq * (double)i / sample_rate);
    }

    unsigned num_bins = N / 2 + 1;
    double *out = malloc(num_bins * sizeof(double));
    fn(sig, N, out);

    unsigned peak_bin = (unsigned)(freq * N / sample_rate);

    unsigned actual_peak = 1;
    for (unsigned k = 2; k < num_bins; k++) {
        if (out[k] > out[actual_peak]) actual_peak = k;
    }

    free(out);
    free(sig);
    return (actual_peak == peak_bin);
}

/** Two sinusoids should produce two peaks. */
static int _test_spectrum_two_sines(spectrum_fn_t fn)
{
    unsigned N = 2048;
    double sample_rate = 2048.0;
    double freq1 = 200.0;
    double freq2 = 500.0;

    double *sig = malloc(N * sizeof(double));
    for (unsigned i = 0; i < N; i++) {
        sig[i] = sin(2.0 * M_PI * freq1 * (double)i / sample_rate)
               + sin(2.0 * M_PI * freq2 * (double)i / sample_rate);
    }

    unsigned num_bins = N / 2 + 1;
    double *out = malloc(num_bins * sizeof(double));
    fn(sig, N, out);

    unsigned bin1 = (unsigned)(freq1 * N / sample_rate);
    unsigned bin2 = (unsigned)(freq2 * N / sample_rate);

    int ok = 1;
    ok &= (out[bin1] > out[bin1 - 2] * 10.0);
    ok &= (out[bin1] > out[bin1 + 2] * 10.0);
    ok &= (out[bin2] > out[bin2 - 2] * 10.0);
    ok &= (out[bin2] > out[bin2 + 2] * 10.0);

    free(out);
    free(sig);
    return ok;
}

/** An all-zeros signal should produce zero output everywhere. */
static int _test_spectrum_zeros(spectrum_fn_t fn)
{
    unsigned N = 512;
    double *sig = calloc(N, sizeof(double));

    unsigned num_bins = N / 2 + 1;
    double *out = malloc(num_bins * sizeof(double));
    fn(sig, N, out);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(out[k], 0.0, 1e-15);
    }

    free(out);
    free(sig);
    return ok;
}

/** Output should be non-negative for any input. */
static int _test_spectrum_non_negative(spectrum_fn_t fn)
{
    unsigned N = 512;
    double *sig = malloc(N * sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        sig[i] = sin(2.0 * M_PI * 7.0 * (double)i / (double)N)
               - 0.5 * cos(2.0 * M_PI * 13.0 * (double)i / (double)N);
    }

    unsigned num_bins = N / 2 + 1;
    double *out = malloc(num_bins * sizeof(double));
    fn(sig, N, out);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= (out[k] >= 0.0);
    }

    free(out);
    free(sig);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_magnitude_spectrum()
 *
 * These tests verify that the magnitude spectrum correctly identifies
 * frequency content in known signals.  The key mathematical properties:
 *
 *   - A pure sine at frequency f concentrates energy in bin k = f*N/sr
 *   - DC offset appears only in bin 0
 *   - Parseval's theorem: sum of |X(k)|^2 / N == sum of x(n)^2
 * -----------------------------------------------------------------------*/

static int test_mag_spectrum_single_sine(void) { return _test_spectrum_single_sine(MD_magnitude_spectrum); }
static int test_mag_spectrum_two_sines(void) { return _test_spectrum_two_sines(MD_magnitude_spectrum); }

/** A DC signal (constant value) should have energy only at bin 0. */
static int test_mag_spectrum_dc(void)
{
    unsigned N = 256;
    double *sig = malloc(N * sizeof(double));
    for (unsigned i = 0; i < N; i++) {
        sig[i] = 3.0;
    }

    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    int ok = 1;
    /* DC bin should have magnitude = N * amplitude = 256 * 3.0 = 768.0 */
    ok &= approx_equal(mag[0], (double)N * 3.0, 1e-10);

    /* All other bins should be zero */
    for (unsigned k = 1; k < num_bins; k++) {
        ok &= approx_equal(mag[k], 0.0, 1e-10);
    }

    free(mag);
    free(sig);

    return ok;
}

static int test_mag_spectrum_zeros(void) { return _test_spectrum_zeros(MD_magnitude_spectrum); }

/**
 * A unit impulse (delta function) has a flat magnitude spectrum.
 * All bins should have |X(k)| = 1.0 because the DFT of delta[n] = 1
 * for all k.
 */
static int test_mag_spectrum_impulse(void)
{
    unsigned N = 256;
    double *sig = calloc(N, sizeof(double));
    sig[0] = 1.0;  /* unit impulse at sample 0 */

    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(mag[k], 1.0, 1e-10);
    }

    free(mag);
    free(sig);

    return ok;
}

/**
 * Parseval's theorem:  sum(|x[n]|^2) == sum(|X[k]|^2) / N
 *
 * For a real signal, the full N-point power sum in the frequency domain
 * equals: |X[0]|^2 + 2*sum(|X[k]|^2 for k=1..N/2-1) + |X[N/2]|^2,
 * all divided by N.
 */
static int test_mag_spectrum_parseval(void)
{
    unsigned N = 1024;
    double *sig = malloc(N * sizeof(double));

    /* Generate a signal with multiple frequency components */
    for (unsigned i = 0; i < N; i++) {
        sig[i] = 0.5 * sin(2.0 * M_PI * 50.0 * (double)i / (double)N)
               + 0.3 * cos(2.0 * M_PI * 120.0 * (double)i / (double)N)
               + 0.2;
    }

    /* Time-domain energy */
    double time_energy = 0.0;
    for (unsigned i = 0; i < N; i++) {
        time_energy += sig[i] * sig[i];
    }

    /* Frequency-domain energy via Parseval's theorem */
    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    double freq_energy = mag[0] * mag[0];           /* DC bin */
    for (unsigned k = 1; k < N / 2; k++) {
        freq_energy += 2.0 * mag[k] * mag[k];      /* interior bins appear twice */
    }
    freq_energy += mag[N / 2] * mag[N / 2];         /* Nyquist bin */
    freq_energy /= (double)N;

    free(mag);
    free(sig);

    return approx_equal(time_energy, freq_energy, 1e-8);
}

/** Calling with different N values should work (tests FFT plan re-caching). */
static int test_mag_spectrum_different_lengths(void)
{
    int ok = 1;

    /* First call with N = 512 */
    {
        unsigned N = 512;
        double *sig = calloc(N, sizeof(double));
        sig[0] = 1.0;

        unsigned num_bins = N / 2 + 1;
        double *mag = malloc(num_bins * sizeof(double));
        MD_magnitude_spectrum(sig, N, mag);

        /* Impulse: all bins should be 1.0 */
        for (unsigned k = 0; k < num_bins; k++) {
            ok &= approx_equal(mag[k], 1.0, 1e-10);
        }

        free(mag);
        free(sig);
    }

    /* Second call with N = 256 (forces plan rebuild) */
    {
        unsigned N = 256;
        double *sig = calloc(N, sizeof(double));
        sig[0] = 1.0;

        unsigned num_bins = N / 2 + 1;
        double *mag = malloc(num_bins * sizeof(double));
        MD_magnitude_spectrum(sig, N, mag);

        for (unsigned k = 0; k < num_bins; k++) {
            ok &= approx_equal(mag[k], 1.0, 1e-10);
        }

        free(mag);
        free(sig);
    }

    return ok;
}

/**
 * The magnitude spectrum should be non-negative for any input.
 * Test with a signal containing both positive and negative values.
 */
static int test_mag_spectrum_non_negative(void) { return _test_spectrum_non_negative(MD_magnitude_spectrum); }

/* -----------------------------------------------------------------------
 * Tests for MD_power_spectral_density()
 *
 * These tests mirror the magnitude spectrum tests but verify power
 * (|X(k)|^2 / N) instead of amplitude (|X(k)|).  Key properties:
 *
 *   - PSD[k] = (magnitude[k])^2 / N
 *   - Parseval: PSD[0] + 2*sum(PSD[1..N/2-1]) + PSD[N/2] = sum(x[n]^2)
 *   - PSD is always >= 0 for any input
 * -----------------------------------------------------------------------*/

static int test_psd_single_sine(void) { return _test_spectrum_single_sine(MD_power_spectral_density); }
static int test_psd_two_sines(void) { return _test_spectrum_two_sines(MD_power_spectral_density); }

/**
 * A DC signal (constant value) concentrates all power in bin 0.
 * For a DC signal x[n] = A, the DFT gives X(0) = N*A and X(k) = 0 for k > 0.
 * So PSD[0] = |N*A|^2 / N = N * A^2.
 */
static int test_psd_dc(void)
{
    unsigned N = 256;
    double amp = 3.0;
    double *sig = malloc(N * sizeof(double));
    for (unsigned i = 0; i < N; i++) {
        sig[i] = amp;
    }

    unsigned num_bins = N / 2 + 1;
    double *psd = malloc(num_bins * sizeof(double));
    MD_power_spectral_density(sig, N, psd);

    int ok = 1;
    /* PSD[0] = (N * A)^2 / N = N * A^2 */
    ok &= approx_equal(psd[0], (double)N * amp * amp, 1e-10);

    /* All other bins should be zero */
    for (unsigned k = 1; k < num_bins; k++) {
        ok &= approx_equal(psd[k], 0.0, 1e-10);
    }

    free(psd);
    free(sig);

    return ok;
}

static int test_psd_zeros(void) { return _test_spectrum_zeros(MD_power_spectral_density); }

/**
 * A unit impulse (delta function) has a flat PSD.
 * The DFT of delta[n] = 1 for all k, so PSD[k] = 1^2 / N = 1/N.
 */
static int test_psd_impulse(void)
{
    unsigned N = 256;
    double *sig = calloc(N, sizeof(double));
    sig[0] = 1.0;

    unsigned num_bins = N / 2 + 1;
    double *psd = malloc(num_bins * sizeof(double));
    MD_power_spectral_density(sig, N, psd);

    int ok = 1;
    double expected = 1.0 / (double)N;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(psd[k], expected, 1e-15);
    }

    free(psd);
    free(sig);

    return ok;
}

/**
 * Parseval's theorem for the periodogram:
 *   PSD[0] + 2 * sum(PSD[1..N/2-1]) + PSD[N/2] = sum(x[n]^2)
 *
 * The one-sided PSD sum equals the total signal energy (not average power).
 * This is because PSD[k] = |X(k)|^2 / N, and Parseval gives
 * sum(|X(k)|^2) = N * sum(x[n]^2), so the one-sided reconstruction
 * recovers the full time-domain energy.
 */
static int test_psd_parseval(void)
{
    unsigned N = 1024;
    double *sig = malloc(N * sizeof(double));

    /* Generate a signal with multiple frequency components */
    for (unsigned i = 0; i < N; i++) {
        sig[i] = 0.5 * sin(2.0 * M_PI * 50.0 * (double)i / (double)N)
               + 0.3 * cos(2.0 * M_PI * 120.0 * (double)i / (double)N)
               + 0.2;
    }

    /* Time-domain energy: sum(x[n]^2) */
    double time_energy = 0.0;
    for (unsigned i = 0; i < N; i++) {
        time_energy += sig[i] * sig[i];
    }

    /* Frequency-domain energy via one-sided PSD sum */
    unsigned num_bins = N / 2 + 1;
    double *psd = malloc(num_bins * sizeof(double));
    MD_power_spectral_density(sig, N, psd);

    double freq_energy = psd[0];                     /* DC bin */
    for (unsigned k = 1; k < N / 2; k++) {
        freq_energy += 2.0 * psd[k];                /* interior bins (doubled) */
    }
    freq_energy += psd[N / 2];                       /* Nyquist bin */

    free(psd);
    free(sig);

    return approx_equal(time_energy, freq_energy, 1e-8);
}

/** Calling with different N values should work (tests FFT plan re-caching). */
static int test_psd_different_lengths(void)
{
    int ok = 1;

    /* First call with N = 512 */
    {
        unsigned N = 512;
        double *sig = calloc(N, sizeof(double));
        sig[0] = 1.0;

        unsigned num_bins = N / 2 + 1;
        double *psd = malloc(num_bins * sizeof(double));
        MD_power_spectral_density(sig, N, psd);

        /* Impulse: PSD[k] = 1/N for all k */
        double expected = 1.0 / (double)N;
        for (unsigned k = 0; k < num_bins; k++) {
            ok &= approx_equal(psd[k], expected, 1e-15);
        }

        free(psd);
        free(sig);
    }

    /* Second call with N = 256 (forces plan rebuild) */
    {
        unsigned N = 256;
        double *sig = calloc(N, sizeof(double));
        sig[0] = 1.0;

        unsigned num_bins = N / 2 + 1;
        double *psd = malloc(num_bins * sizeof(double));
        MD_power_spectral_density(sig, N, psd);

        double expected = 1.0 / (double)N;
        for (unsigned k = 0; k < num_bins; k++) {
            ok &= approx_equal(psd[k], expected, 1e-15);
        }

        free(psd);
        free(sig);
    }

    return ok;
}

static int test_psd_non_negative(void) { return _test_spectrum_non_negative(MD_power_spectral_density); }

/* -----------------------------------------------------------------------
 * Tests for MD_phase_spectrum()
 *
 * Phase is well-defined only at bins with significant magnitude.  All
 * tests use exact-integer-bin frequencies (sample_rate == N) so FFTW
 * produces bit-exact complex values with no spectral leakage.
 * -----------------------------------------------------------------------*/

/**
 * A pure cosine at an integer bin k0 has zero phase at that bin.
 *
 * cos(2pi*k0*n/N) = Re[exp(j*2pi*k0*n/N)], so X(k0) is real and
 * positive, giving atan2(0, positive) == 0.
 */
static int test_phase_spectrum_cosine(void)
{
    unsigned N = 1024;
    unsigned k0 = 100;
    double *sig = malloc(N * sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        sig[i] = cos(2.0 * M_PI * (double)k0 * (double)i / (double)N);
    }

    unsigned num_bins = N / 2 + 1;
    double *phase = malloc(num_bins * sizeof(double));
    MD_phase_spectrum(sig, N, phase);

    int ok = approx_equal(phase[k0], 0.0, 1e-9);

    free(phase);
    free(sig);
    return ok;
}

/**
 * A pure sine at an integer bin k0 has phase -pi/2 at that bin.
 *
 * sin(2pi*k0*n/N) = -Im[exp(j*2pi*k0*n/N)], so X(k0) is purely
 * negative-imaginary, giving atan2(-N/2, 0) == -pi/2.
 */
static int test_phase_spectrum_sine(void)
{
    unsigned N = 1024;
    unsigned k0 = 100;
    double *sig = malloc(N * sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        sig[i] = sin(2.0 * M_PI * (double)k0 * (double)i / (double)N);
    }

    unsigned num_bins = N / 2 + 1;
    double *phase = malloc(num_bins * sizeof(double));
    MD_phase_spectrum(sig, N, phase);

    int ok = approx_equal(phase[k0], -M_PI / 2.0, 1e-9);

    free(phase);
    free(sig);
    return ok;
}

/**
 * A unit impulse at n=0 has zero phase everywhere.
 *
 * The DFT of delta[n] is X(k) = 1 for all k (real, positive),
 * so atan2(0, 1) = 0 for every bin.
 */
static int test_phase_spectrum_impulse_at_zero(void)
{
    unsigned N = 512;
    double *sig = calloc(N, sizeof(double));
    sig[0] = 1.0;

    unsigned num_bins = N / 2 + 1;
    double *phase = malloc(num_bins * sizeof(double));
    MD_phase_spectrum(sig, N, phase);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(phase[k], 0.0, 1e-12);
    }

    free(phase);
    free(sig);
    return ok;
}

/**
 * A unit impulse at n=1 produces linear phase: phi(k) = -2*pi*k/N.
 *
 * By the DFT shift theorem, delaying a signal by d samples multiplies
 * each frequency bin by exp(-j*2*pi*k*d/N), adding a linear phase
 * ramp of slope -2*pi*d/N radians per bin.
 */
static int test_phase_spectrum_impulse_at_one(void)
{
    unsigned N = 512;
    double *sig = calloc(N, sizeof(double));
    sig[1] = 1.0;   /* unit impulse at n = 1 (delay d = 1) */

    unsigned num_bins = N / 2 + 1;
    double *phase = malloc(num_bins * sizeof(double));
    MD_phase_spectrum(sig, N, phase);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        double expected = -2.0 * M_PI * (double)k / (double)N;
        /* Use a wrap-aware comparison: -pi and +pi are the same angle.
         * At the Nyquist bin (k = N/2), expected = -pi.  FFTW returns +pi
         * because the tiny floating-point imaginary part can have either
         * sign, so carg() may land on either boundary of [-pi, pi]. */
        double diff = fabs(phase[k] - expected);
        if (diff > M_PI) diff = 2.0 * M_PI - diff;  /* wrap to [0, pi] */
        ok &= (diff < 1e-12);
    }

    free(phase);
    free(sig);
    return ok;
}

/**
 * An all-zeros signal should produce zero phase everywhere.
 *
 * carg(0 + 0j) == 0 by convention (atan2(0, 0) == 0 in C).
 */
static int test_phase_spectrum_zeros(void)
{
    unsigned N = 256;
    double *sig = calloc(N, sizeof(double));

    unsigned num_bins = N / 2 + 1;
    double *phase = malloc(num_bins * sizeof(double));
    MD_phase_spectrum(sig, N, phase);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(phase[k], 0.0, 1e-15);
    }

    free(phase);
    free(sig);
    return ok;
}

/**
 * Calling MD_phase_spectrum() with two different N values tests
 * that the FFT plan cache is correctly rebuilt when N changes.
 */
static int test_phase_spectrum_different_lengths(void)
{
    int ok = 1;

    /* First call: cosine at bin 50 in N=512 */
    {
        unsigned N = 512;
        unsigned k0 = 50;
        double *sig = malloc(N * sizeof(double));
        for (unsigned i = 0; i < N; i++) {
            sig[i] = cos(2.0 * M_PI * (double)k0 * (double)i / (double)N);
        }
        unsigned num_bins = N / 2 + 1;
        double *phase = malloc(num_bins * sizeof(double));
        MD_phase_spectrum(sig, N, phase);
        ok &= approx_equal(phase[k0], 0.0, 1e-9);
        free(phase);
        free(sig);
    }

    /* Second call: cosine at bin 200 in N=2048 (forces plan rebuild) */
    {
        unsigned N = 2048;
        unsigned k0 = 200;
        double *sig = malloc(N * sizeof(double));
        for (unsigned i = 0; i < N; i++) {
            sig[i] = cos(2.0 * M_PI * (double)k0 * (double)i / (double)N);
        }
        unsigned num_bins = N / 2 + 1;
        double *phase = malloc(num_bins * sizeof(double));
        MD_phase_spectrum(sig, N, phase);
        ok &= approx_equal(phase[k0], 0.0, 1e-9);
        free(phase);
        free(sig);
    }

    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_gcc() and MD_get_delay() -- GCC-PHAT
 *
 * These tests create known delayed signals and verify the algorithm
 * correctly recovers the delay.
 * -----------------------------------------------------------------------*/

/** Detect a positive delay using PHAT weighting. */
static int test_gcc_phat_positive_delay(void)
{
    int true_delay = 5;
    unsigned N = 4096;
    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));

    /* Generate a test signal */
    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }
    delay_signal(siga, sigb, N, true_delay);

    int delay = MD_get_delay(siga, sigb, N, nullptr, 50, PHAT);
    free(sigb);
    free(siga);

    return (delay == true_delay);
}

/** Detect a negative delay using PHAT weighting. */
static int test_gcc_phat_negative_delay(void)
{
    int true_delay = -7;
    unsigned N = 8192;
    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        siga[i] = 10.0 * sin(2.0 * M_PI * (double)i / 128.0);
    }
    delay_signal(siga, sigb, N, true_delay);

    int delay = MD_get_delay(siga, sigb, N, nullptr, 50, PHAT);
    free(sigb);
    free(siga);

    return (delay == true_delay);
}

/** Zero delay should give zero. */
static int test_gcc_phat_zero_delay(void)
{
    unsigned N = 4096;
    double *siga = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }

    /* Both signals are identical -- delay should be 0 */
    int delay = MD_get_delay(siga, siga, N, nullptr, 50, PHAT);
    free(siga);

    return (delay == 0);
}

/** Test with SIMP (non-PHAT) weighting. */
static int test_gcc_simp_delay(void)
{
    int true_delay = 3;
    unsigned N = 4096;
    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }
    delay_signal(siga, sigb, N, true_delay);

    int delay = MD_get_delay(siga, sigb, N, nullptr, 50, SIMP);
    free(sigb);
    free(siga);

    return (delay == true_delay);
}

/** MD_get_delay should return entropy when ent is not null. */
static int test_gcc_returns_entropy(void)
{
    unsigned N = 4096;
    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }
    delay_signal(siga, sigb, N, 3);

    double ent = -1.0;
    MD_get_delay(siga, sigb, N, &ent, 50, PHAT);
    free(sigb);
    free(siga);

    /* Entropy should be a valid value between 0 and 1 */
    return (ent >= 0.0 && ent <= 1.0);
}

/** Test MD_get_multiple_delays agrees with manual windowed MD_get_delay.
 *
 * MD_get_multiple_delays is a convenience wrapper that:
 *   1. Applies a Hanning window to each signal.
 *   2. Calls MD_get_delay on each windowed pair.
 *
 * This test verifies that its output matches what we get by doing
 * those same steps manually.
 */
static int test_gcc_multiple_delays(void)
{
    unsigned N = 8192;
    int delays[] = {-7, 5, -2};

    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));
    double *sigc = calloc(N, sizeof(double));
    double *sigd = calloc(N, sizeof(double));

    /* Generate reference signal */
    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }
    delay_signal(siga, sigb, N, delays[0]);
    delay_signal(siga, sigc, N, delays[1]);
    delay_signal(siga, sigd, N, delays[2]);

    /* Get results from MD_get_multiple_delays */
    const double *sigs[] = {siga, sigb, sigc, sigd};
    int multi_results[3];
    MD_get_multiple_delays(sigs, 4, N, 50, PHAT, multi_results);

    /* Manually window and call MD_get_delay for comparison */
    double *hann = calloc(N, sizeof(double));
    MD_Gen_Hann_Win(hann, N);

    double *w_ref = calloc(N, sizeof(double));
    double *w_sig = calloc(N, sizeof(double));

    int manual_results[3];
    const double *other[] = {sigb, sigc, sigd};
    for (int k = 0; k < 3; k++) {
        for (unsigned i = 0; i < N; i++) {
            w_ref[i] = siga[i] * hann[i];
            w_sig[i] = other[k][i] * hann[i];
        }
        manual_results[k] = MD_get_delay(w_ref, w_sig, N, nullptr, 50, PHAT);
    }

    int ok = 1;
    ok &= (multi_results[0] == manual_results[0]);
    ok &= (multi_results[1] == manual_results[1]);
    ok &= (multi_results[2] == manual_results[2]);

    free(w_sig);
    free(w_ref);
    free(hann);
    free(sigd);
    free(sigc);
    free(sigb);
    free(siga);

    return ok;
}

/** The gcc output array should have its peak at the right place. */
static int test_gcc_lagvals_structure(void)
{
    unsigned N = 4096;
    double *siga = calloc(N, sizeof(double));
    double *lagvals = calloc(N, sizeof(double));

    for (unsigned i = 0; i < N; i++) {
        siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
    }

    /* Auto-correlation: the peak should be at the zero-lag position */
    MD_gcc(siga, siga, N, lagvals, PHAT);

    unsigned center = (unsigned)ceil((double)N / 2.0);

    /* Find the index of the maximum */
    unsigned max_i = 0;
    double max_v = lagvals[0];
    for (unsigned i = 1; i < N; i++) {
        if (lagvals[i] > max_v) {
            max_v = lagvals[i];
            max_i = i;
        }
    }

    free(lagvals);
    free(siga);

    /* The peak should be at the center (zero-lag) */
    return (max_i == center);
}

/** Verify delay detection works after changing signal length (tests FFT plan caching). */
static int test_gcc_different_lengths(void)
{
    int ok = 1;

    /* First call with N=4096 */
    {
        int true_delay = 3;
        unsigned N = 4096;
        double *siga = calloc(N, sizeof(double));
        double *sigb = calloc(N, sizeof(double));
        for (unsigned i = 0; i < N; i++)
            siga[i] = sin(2.0 * M_PI * (double)i / 64.0);
        delay_signal(siga, sigb, N, true_delay);
        int d = MD_get_delay(siga, sigb, N, nullptr, 50, PHAT);
        ok &= (d == true_delay);
        free(sigb);
        free(siga);
    }

    /* Second call with N=2048 (forces FFT plan rebuild) */
    {
        int true_delay = -4;
        unsigned N = 2048;
        double *siga = calloc(N, sizeof(double));
        double *sigb = calloc(N, sizeof(double));
        for (unsigned i = 0; i < N; i++)
            siga[i] = sin(2.0 * M_PI * (double)i / 32.0);
        delay_signal(siga, sigb, N, true_delay);
        int d = MD_get_delay(siga, sigb, N, nullptr, 50, PHAT);
        ok &= (d == true_delay);
        free(sigb);
        free(siga);
    }

    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for BiQuad filters
 *
 * We verify each filter type by feeding a known signal through it and
 * checking that the expected frequency components are preserved or
 * attenuated.
 * -----------------------------------------------------------------------*/

/**
 * Helper: measure the RMS amplitude of a sine wave after filtering.
 *
 * Generates a sine wave at the given frequency, passes it through
 * the biquad filter, and returns the RMS of the steady-state output
 * (skipping the initial transient).
 */
static double measure_filter_response(biquad *b, double freq, double srate,
                                      unsigned total_samples, unsigned skip)
{
    double sum_sq = 0.0;
    unsigned count = 0;

    for (unsigned i = 0; i < total_samples; i++) {
        double sample = sin(2.0 * M_PI * freq * (double)i / srate);
        double out = BiQuad(sample, b);
        if (i >= skip) {
            sum_sq += out * out;
            count++;
        }
    }

    return sqrt(sum_sq / (double)count);
}

/** BiQuad_new should return nullptr for an invalid filter type. */
static int test_biquad_invalid_type(void)
{
    biquad *b = BiQuad_new(999, 0.0, 1000.0, 44100.0, 1.0);
    return (b == nullptr);
}

/** Low-pass filter: should pass low frequencies, attenuate high. */
static int test_biquad_lpf(void)
{
    double srate = 44100.0;
    double cutoff = 1000.0;
    unsigned nsamples = 44100;  /* 1 second */
    unsigned skip = 4410;       /* skip first 0.1s transient */

    /* Test with a frequency well below the cutoff */
    biquad *b_low = BiQuad_new(LPF, 0.0, cutoff, srate, 1.0);
    double rms_low = measure_filter_response(b_low, 200.0, srate, nsamples, skip);
    free(b_low);

    /* Test with a frequency well above the cutoff */
    biquad *b_high = BiQuad_new(LPF, 0.0, cutoff, srate, 1.0);
    double rms_high = measure_filter_response(b_high, 10000.0, srate, nsamples, skip);
    free(b_high);

    /* The low frequency should pass through much louder than the high frequency.
     * The input sine has RMS = 1/sqrt(2) ≈ 0.707.
     * A good LPF should keep the low freq near that and attenuate the high freq. */
    return (rms_low > 0.5 && rms_high < 0.1);
}

/** High-pass filter: should pass high frequencies, attenuate low. */
static int test_biquad_hpf(void)
{
    double srate = 44100.0;
    double cutoff = 1000.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    biquad *b_low = BiQuad_new(HPF, 0.0, cutoff, srate, 1.0);
    double rms_low = measure_filter_response(b_low, 100.0, srate, nsamples, skip);
    free(b_low);

    biquad *b_high = BiQuad_new(HPF, 0.0, cutoff, srate, 1.0);
    double rms_high = measure_filter_response(b_high, 10000.0, srate, nsamples, skip);
    free(b_high);

    return (rms_high > 0.5 && rms_low < 0.1);
}

/** Band-pass filter: should pass the centre frequency, attenuate others. */
static int test_biquad_bpf(void)
{
    double srate = 44100.0;
    double center = 2000.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    /* At centre frequency */
    biquad *b_center = BiQuad_new(BPF, 0.0, center, srate, 1.0);
    double rms_center = measure_filter_response(b_center, center, srate, nsamples, skip);
    free(b_center);

    /* Far below */
    biquad *b_low = BiQuad_new(BPF, 0.0, center, srate, 1.0);
    double rms_low = measure_filter_response(b_low, 100.0, srate, nsamples, skip);
    free(b_low);

    /* Far above */
    biquad *b_high = BiQuad_new(BPF, 0.0, center, srate, 1.0);
    double rms_high = measure_filter_response(b_high, 15000.0, srate, nsamples, skip);
    free(b_high);

    return (rms_center > rms_low && rms_center > rms_high && rms_center > 0.3);
}

/** Notch filter: should remove the centre frequency, pass others. */
static int test_biquad_notch(void)
{
    double srate = 44100.0;
    double notch_freq = 2000.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    /* At the notch frequency (should be attenuated) */
    biquad *b_notch = BiQuad_new(NOTCH, 0.0, notch_freq, srate, 1.0);
    double rms_notch = measure_filter_response(b_notch, notch_freq, srate, nsamples, skip);
    free(b_notch);

    /* Away from the notch (should pass through) */
    biquad *b_pass = BiQuad_new(NOTCH, 0.0, notch_freq, srate, 1.0);
    double rms_pass = measure_filter_response(b_pass, 500.0, srate, nsamples, skip);
    free(b_pass);

    return (rms_notch < 0.1 && rms_pass > 0.5);
}

/** Peaking EQ with positive gain should boost the centre frequency. */
static int test_biquad_peq_boost(void)
{
    double srate = 44100.0;
    double center = 2000.0;
    double gain_db = 12.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    biquad *b = BiQuad_new(PEQ, gain_db, center, srate, 1.0);
    double rms_boosted = measure_filter_response(b, center, srate, nsamples, skip);
    free(b);

    /* Input RMS is about 0.707; with +12dB boost, output should be much higher */
    return (rms_boosted > 1.5);
}

/** Low shelf with positive gain should boost low frequencies. */
static int test_biquad_low_shelf(void)
{
    double srate = 44100.0;
    double shelf_freq = 1000.0;
    double gain_db = 12.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    /* Low frequency should be boosted */
    biquad *b_low = BiQuad_new(LSH, gain_db, shelf_freq, srate, 1.0);
    double rms_low = measure_filter_response(b_low, 100.0, srate, nsamples, skip);
    free(b_low);

    /* High frequency should be relatively unaffected */
    biquad *b_high = BiQuad_new(LSH, gain_db, shelf_freq, srate, 1.0);
    double rms_high = measure_filter_response(b_high, 10000.0, srate, nsamples, skip);
    free(b_high);

    return (rms_low > rms_high && rms_low > 1.0);
}

/** High shelf with positive gain should boost high frequencies. */
static int test_biquad_high_shelf(void)
{
    double srate = 44100.0;
    double shelf_freq = 1000.0;
    double gain_db = 12.0;
    unsigned nsamples = 44100;
    unsigned skip = 4410;

    /* High frequency should be boosted */
    biquad *b_high = BiQuad_new(HSH, gain_db, shelf_freq, srate, 1.0);
    double rms_high = measure_filter_response(b_high, 10000.0, srate, nsamples, skip);
    free(b_high);

    /* Low frequency should be relatively unaffected */
    biquad *b_low = BiQuad_new(HSH, gain_db, shelf_freq, srate, 1.0);
    double rms_low = measure_filter_response(b_low, 100.0, srate, nsamples, skip);
    free(b_low);

    return (rms_high > rms_low && rms_high > 1.0);
}

/** Filtering a DC signal through a HPF should output zero (steady state). */
static int test_biquad_hpf_dc_rejection(void)
{
    biquad *b = BiQuad_new(HPF, 0.0, 100.0, 44100.0, 1.0);
    double last_output = 0.0;

    /* Feed a constant (DC) signal through the high-pass filter */
    for (unsigned i = 0; i < 44100; i++) {
        last_output = BiQuad(1.0, b);
    }
    free(b);

    /* After settling, output should be essentially zero */
    return approx_equal(last_output, 0.0, 1e-6);
}

/* -----------------------------------------------------------------------
 * Tests for MD_sine_wave()
 * -----------------------------------------------------------------------*/

/** output[0] must be 0 for any amplitude/freq because sin(0) = 0. */
static int test_sine_starts_at_zero(void)
{
    double out[64];
    MD_sine_wave(out, 64, 3.7, 440.0, 44100.0);
    return approx_equal(out[0], 0.0, 1e-12);
}

/** At i = sample_rate / (4 * freq) the sine reaches its peak (sin(π/2) = 1). */
static int test_sine_quarter_period(void)
{
    double sample_rate = 16000.0;
    double freq = 1000.0;
    double amplitude = 2.5;
    unsigned N = (unsigned)(sample_rate / freq) * 4;
    double out[N];
    MD_sine_wave(out, N, amplitude, freq, sample_rate);
    unsigned quarter = (unsigned)(sample_rate / (4.0 * freq));
    return approx_equal(out[quarter], amplitude, 1e-9);
}

/** At i = sample_rate / freq the sine completes one full period (sin(2π) ≈ 0). */
static int test_sine_full_period(void)
{
    double sample_rate = 16000.0;
    double freq = 500.0;
    double amplitude = 1.0;
    unsigned period = (unsigned)(sample_rate / freq);
    double out[period + 1];
    MD_sine_wave(out, period + 1, amplitude, freq, sample_rate);
    return approx_equal(out[period], 0.0, 1e-9);
}

/** Magnitude spectrum of a 1 kHz sine (N=1024, fs=16 kHz) peaks at bin k=64. */
static int test_sine_spectrum_peak(void)
{
    unsigned N = 1024;
    double sample_rate = 16000.0;
    double freq = 1000.0;                         /* bin 64: 64 * 16000/1024 */
    unsigned expected_bin = (unsigned)(freq * N / sample_rate);  /* = 64 */

    double sig[N];
    MD_sine_wave(sig, N, 1.0, freq, sample_rate);

    unsigned num_bins = N / 2 + 1;
    double mag[num_bins];
    MD_magnitude_spectrum(sig, N, mag);

    unsigned peak_bin = 0;
    for (unsigned k = 1; k < num_bins; k++) {
        if (mag[k] > mag[peak_bin])
            peak_bin = k;
    }
    return (peak_bin == expected_bin);
}

/** Negative amplitude: output[0]=0, quarter-period value ≈ amplitude. */
static int test_sine_amplitude_negative(void)
{
    double sample_rate = 44100.0;
    double freq = 441.0;
    double amplitude = -2.0;
    unsigned N = 1024;
    double out[N];
    MD_sine_wave(out, N, amplitude, freq, sample_rate);
    unsigned quarter = (unsigned)(sample_rate / (4.0 * freq));
    int ok = approx_equal(out[0], 0.0, 1e-12);
    ok &= approx_equal(out[quarter], amplitude, 1e-6);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_white_noise()
 * -----------------------------------------------------------------------*/

/** Sample mean of white noise should be near zero. */
static int test_white_noise_mean_near_zero(void)
{
    unsigned N = 100000;
    double *noise = malloc(N * sizeof(double));
    MD_white_noise(noise, N, 1.0, 42);

    double sum = 0.0;
    for (unsigned i = 0; i < N; i++) sum += noise[i];
    double mean = sum / (double)N;

    free(noise);
    return (fabs(mean) < 0.02);
}

/** Sample standard deviation should match the requested amplitude. */
static int test_white_noise_stddev_matches_amplitude(void)
{
    unsigned N = 100000;
    double amplitude = 2.5;
    double *noise = malloc(N * sizeof(double));
    MD_white_noise(noise, N, amplitude, 123);

    double sum = 0.0, sum_sq = 0.0;
    for (unsigned i = 0; i < N; i++) {
        sum += noise[i];
        sum_sq += noise[i] * noise[i];
    }
    double mean = sum / (double)N;
    double variance = sum_sq / (double)N - mean * mean;
    double stddev = sqrt(variance);

    free(noise);
    return approx_equal(stddev, amplitude, 0.05);
}

/** Same seed must produce identical output. */
static int test_white_noise_reproducible(void)
{
    unsigned N = 1024;
    double *a = malloc(N * sizeof(double));
    double *b = malloc(N * sizeof(double));
    MD_white_noise(a, N, 1.0, 99);
    MD_white_noise(b, N, 1.0, 99);

    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        ok &= (a[i] == b[i]);
    }
    free(b);
    free(a);
    return ok;
}

/** Different seeds must produce different output. */
static int test_white_noise_different_seeds(void)
{
    unsigned N = 1024;
    double *a = malloc(N * sizeof(double));
    double *b = malloc(N * sizeof(double));
    MD_white_noise(a, N, 1.0, 1);
    MD_white_noise(b, N, 1.0, 2);

    int differs = 0;
    for (unsigned i = 0; i < N; i++) {
        if (a[i] != b[i]) { differs = 1; break; }
    }
    free(b);
    free(a);
    return differs;
}

/** PSD of white noise should be roughly flat (no bin > 10x the mean). */
static int test_white_noise_flat_spectrum(void)
{
    unsigned N = 4096;
    double *noise = malloc(N * sizeof(double));
    MD_white_noise(noise, N, 1.0, 77);

    unsigned num_bins = N / 2 + 1;
    double *psd = malloc(num_bins * sizeof(double));
    MD_power_spectral_density(noise, N, psd);

    /* Compute mean PSD (skip DC bin 0) */
    double sum = 0.0;
    for (unsigned k = 1; k < num_bins; k++) sum += psd[k];
    double mean_psd = sum / (double)(num_bins - 1);

    /* No bin should be more than 10x the mean */
    int ok = (mean_psd > 0.0);
    for (unsigned k = 1; k < num_bins; k++) {
        if (psd[k] > mean_psd * 10.0) { ok = 0; break; }
    }

    free(psd);
    free(noise);
    return ok;
}

/** Odd lengths (including N=1) must produce finite values. */
static int test_white_noise_odd_length(void)
{
    double out1[1];
    MD_white_noise(out1, 1, 1.0, 0);
    if (isnan(out1[0]) || isinf(out1[0])) return 0;

    double out7[7];
    MD_white_noise(out7, 7, 1.0, 0);
    for (unsigned i = 0; i < 7; i++) {
        if (isnan(out7[i]) || isinf(out7[i])) return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * Tests for MD_impulse()
 * -----------------------------------------------------------------------*/

/** Unit impulse at position 0: output[0]=1.0, all others zero. */
static int test_impulse_at_zero(void)
{
    unsigned N = 64;
    double out[64];
    MD_impulse(out, N, 1.0, 0);
    if (!approx_equal(out[0], 1.0, 1e-15)) return 0;
    for (unsigned i = 1; i < N; i++) {
        if (!approx_equal(out[i], 0.0, 1e-15)) return 0;
    }
    return 1;
}

/** Impulse at the middle of the buffer. */
static int test_impulse_at_middle(void)
{
    unsigned N = 128;
    double out[128];
    unsigned pos = N / 2;
    MD_impulse(out, N, 1.0, pos);
    for (unsigned i = 0; i < N; i++) {
        double expected = (i == pos) ? 1.0 : 0.0;
        if (!approx_equal(out[i], expected, 1e-15)) return 0;
    }
    return 1;
}

/** Impulse at the last sample (position = N-1). */
static int test_impulse_at_last(void)
{
    unsigned N = 256;
    double out[256];
    MD_impulse(out, N, 1.0, N - 1);
    for (unsigned i = 0; i < N - 1; i++) {
        if (!approx_equal(out[i], 0.0, 1e-15)) return 0;
    }
    return approx_equal(out[N - 1], 1.0, 1e-15);
}

/** Custom amplitude: the spike value should match the given amplitude. */
static int test_impulse_amplitude(void)
{
    unsigned N = 32;
    double out[32];
    MD_impulse(out, N, 5.5, 10);
    int ok = approx_equal(out[10], 5.5, 1e-15);
    ok &= approx_equal(out[9], 0.0, 1e-15);
    ok &= approx_equal(out[11], 0.0, 1e-15);
    return ok;
}

/** Negative amplitude should work (just a negative spike). */
static int test_impulse_amplitude_negative(void)
{
    unsigned N = 64;
    double out[64];
    MD_impulse(out, N, -3.0, 0);
    return approx_equal(out[0], -3.0, 1e-15);
}

/** N=1: a single-sample buffer should contain only the impulse. */
static int test_impulse_single_sample(void)
{
    double out[1];
    MD_impulse(out, 1, 2.0, 0);
    return approx_equal(out[0], 2.0, 1e-15);
}

/** Unit impulse at position 0 has a flat magnitude spectrum. */
static int test_impulse_flat_spectrum(void)
{
    unsigned N = 512;
    double *sig = malloc(N * sizeof(double));
    MD_impulse(sig, N, 1.0, 0);

    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++)
        ok &= approx_equal(mag[k], 1.0, 1e-10);

    free(mag);
    free(sig);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_chirp_linear()
 * -----------------------------------------------------------------------*/

/** output[0] must be 0 for any chirp because sin(0) = 0. */
static int test_chirp_linear_starts_at_zero(void)
{
    double out[1024];
    MD_chirp_linear(out, 1024, 2.5, 200.0, 4000.0, 16000.0);
    return approx_equal(out[0], 0.0, 1e-12);
}

/** The magnitude spectrum of a linear chirp should have energy spread
 *  across the swept frequency range, not concentrated in one bin. */
static int test_chirp_linear_spectrum_spread(void)
{
    unsigned N = 4096;
    double sample_rate = 16000.0;
    double f_start = 1000.0;
    double f_end   = 3000.0;

    double *sig = malloc(N * sizeof(double));
    MD_chirp_linear(sig, N, 1.0, f_start, f_end, sample_rate);

    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    unsigned bin_start = (unsigned)(f_start * N / sample_rate);
    unsigned bin_end   = (unsigned)(f_end   * N / sample_rate);

    /* Find the peak magnitude in the sweep range */
    double peak = 0.0;
    for (unsigned k = bin_start; k <= bin_end; k++)
        if (mag[k] > peak) peak = mag[k];

    /* Count bins with energy above 10% of peak */
    double threshold = peak * 0.1;
    unsigned active_bins = 0;
    for (unsigned k = bin_start; k <= bin_end; k++)
        if (mag[k] > threshold) active_bins++;

    free(mag);
    free(sig);

    /* A chirp should have energy in many bins, not just one. */
    return (active_bins > (bin_end - bin_start) / 2);
}

/** When f_start == f_end, the linear chirp degenerates to a constant-
 *  frequency sine wave and should match MD_sine_wave output. */
static int test_chirp_linear_constant_freq(void)
{
    unsigned N = 1024;
    double sample_rate = 1024.0;
    double freq = 100.0;

    double *chirp = malloc(N * sizeof(double));
    double *sine  = malloc(N * sizeof(double));
    MD_chirp_linear(chirp, N, 1.0, freq, freq, sample_rate);
    MD_sine_wave(sine, N, 1.0, freq, sample_rate);

    int ok = 1;
    for (unsigned i = 0; i < N; i++)
        ok &= approx_equal(chirp[i], sine[i], 1e-10);

    free(sine);
    free(chirp);
    return ok;
}

/** The signal should be bounded by [-amplitude, +amplitude]. */
static int test_chirp_linear_amplitude_bound(void)
{
    unsigned N = 8192;
    double amplitude = 3.7;
    double out[N];
    MD_chirp_linear(out, N, amplitude, 100.0, 5000.0, 44100.0);

    for (unsigned i = 0; i < N; i++)
        if (out[i] > amplitude + 1e-12 || out[i] < -amplitude - 1e-12)
            return 0;
    return 1;
}

/* -----------------------------------------------------------------------
 * Tests for MD_chirp_log()
 * -----------------------------------------------------------------------*/

/** output[0] must be 0 because sin(0) = 0. */
static int test_chirp_log_starts_at_zero(void)
{
    double out[1024];
    MD_chirp_log(out, 1024, 1.0, 20.0, 20000.0, 44100.0);
    return approx_equal(out[0], 0.0, 1e-12);
}

/** The log chirp spectrum should have energy spread across the range. */
static int test_chirp_log_spectrum_spread(void)
{
    unsigned N = 4096;
    double sample_rate = 16000.0;
    double f_start = 500.0;
    double f_end   = 4000.0;

    double *sig = malloc(N * sizeof(double));
    MD_chirp_log(sig, N, 1.0, f_start, f_end, sample_rate);

    unsigned num_bins = N / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    MD_magnitude_spectrum(sig, N, mag);

    unsigned bin_start = (unsigned)(f_start * N / sample_rate);
    unsigned bin_end   = (unsigned)(f_end   * N / sample_rate);

    double peak = 0.0;
    for (unsigned k = bin_start; k <= bin_end; k++)
        if (mag[k] > peak) peak = mag[k];

    double threshold = peak * 0.1;
    unsigned active_bins = 0;
    for (unsigned k = bin_start; k <= bin_end; k++)
        if (mag[k] > threshold) active_bins++;

    free(mag);
    free(sig);

    return (active_bins > (bin_end - bin_start) / 2);
}

/** The signal should be bounded by [-amplitude, +amplitude]. */
static int test_chirp_log_amplitude_bound(void)
{
    unsigned N = 8192;
    double amplitude = 2.0;
    double out[N];
    MD_chirp_log(out, N, amplitude, 100.0, 10000.0, 44100.0);

    for (unsigned i = 0; i < N; i++)
        if (out[i] > amplitude + 1e-12 || out[i] < -amplitude - 1e-12)
            return 0;
    return 1;
}

/* -----------------------------------------------------------------------
 * Tests for MD_square_wave()
 * -----------------------------------------------------------------------*/

/** First half of the period is +amplitude, second half is −amplitude. */
static int test_square_high_low(void)
{
    double sample_rate = 16000.0;
    double freq = 1000.0;
    double amplitude = 2.5;
    /* One full period = 16 samples */
    unsigned period = (unsigned)(sample_rate / freq);
    unsigned N = period;
    double out[N];
    MD_square_wave(out, N, amplitude, freq, sample_rate);

    int ok = 1;
    /* Sample 0: phase = 0 → zero crossing */
    ok &= approx_equal(out[0], 0.0, 1e-12);
    /* Samples 1..7 (first half, excluding zero crossing): +amplitude */
    for (unsigned i = 1; i < period / 2; i++)
        ok &= approx_equal(out[i], amplitude, 1e-12);
    /* Sample 8: phase = π → zero crossing */
    ok &= approx_equal(out[period / 2], 0.0, 1e-12);
    /* Samples 9..15 (second half): −amplitude */
    for (unsigned i = period / 2 + 1; i < period; i++)
        ok &= approx_equal(out[i], -amplitude, 1e-12);
    return ok;
}

/** Spectrum of a square wave: peaks at odd harmonics, near-zero at even. */
static int test_square_spectrum_harmonics(void)
{
    unsigned N = 4096;
    double sample_rate = 16000.0;
    double freq = 250.0;  /* bin 64: 64 * 16000/4096 */

    double sig[N];
    MD_square_wave(sig, N, 1.0, freq, sample_rate);

    unsigned num_bins = N / 2 + 1;
    double mag[num_bins];
    MD_magnitude_spectrum(sig, N, mag);

    unsigned fundamental_bin = (unsigned)(freq * N / sample_rate);
    int ok = 1;
    /* Odd harmonics (1f, 3f, 5f) should have significant energy */
    ok &= (mag[fundamental_bin * 1] > mag[fundamental_bin * 2] * 5.0);
    ok &= (mag[fundamental_bin * 3] > mag[fundamental_bin * 2] * 2.0);
    ok &= (mag[fundamental_bin * 5] > mag[fundamental_bin * 4] * 2.0);
    /* Even harmonics (2f, 4f) should be near zero relative to fundamental */
    ok &= (mag[fundamental_bin * 2] < mag[fundamental_bin] * 0.05);
    ok &= (mag[fundamental_bin * 4] < mag[fundamental_bin] * 0.05);
    return ok;
}

/** Negative amplitude inverts the square wave. */
static int test_square_amplitude_negative(void)
{
    double sample_rate = 16000.0;
    double freq = 1000.0;
    double amplitude = -3.0;
    unsigned period = (unsigned)(sample_rate / freq);
    unsigned N = period;
    double out[N];
    MD_square_wave(out, N, amplitude, freq, sample_rate);

    int ok = 1;
    ok &= approx_equal(out[0], 0.0, 1e-12);
    /* With negative amplitude: first half is negative, second half is positive */
    ok &= approx_equal(out[1], amplitude, 1e-12);
    ok &= approx_equal(out[period / 2 + 1], -amplitude, 1e-12);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_sawtooth_wave()
 * -----------------------------------------------------------------------*/

/** Middle of the period is near 0. */
static int test_sawtooth_midpoint(void)
{
    double sample_rate = 16000.0;
    double freq = 1000.0;
    double amplitude = 2.0;
    unsigned period = (unsigned)(sample_rate / freq);
    unsigned N = period;
    double out[N];
    MD_sawtooth_wave(out, N, amplitude, freq, sample_rate);

    /* Start of period: −amplitude */
    int ok = approx_equal(out[0], -amplitude, 1e-12);
    /* Middle of period: near 0 */
    ok &= approx_equal(out[period / 2], 0.0, 1e-9);
    return ok;
}

/** Spectrum: all harmonics present, decaying approximately as 1/k. */
static int test_sawtooth_spectrum_harmonics(void)
{
    unsigned N = 4096;
    double sample_rate = 16000.0;
    double freq = 250.0;  /* bin 64 */

    double sig[N];
    MD_sawtooth_wave(sig, N, 1.0, freq, sample_rate);

    unsigned num_bins = N / 2 + 1;
    double mag[num_bins];
    MD_magnitude_spectrum(sig, N, mag);

    unsigned fb = (unsigned)(freq * N / sample_rate);
    int ok = 1;
    /* All harmonics 1f..5f should have significant energy */
    for (unsigned k = 1; k <= 5; k++)
        ok &= (mag[fb * k] > mag[fb] * 0.05);
    /* Harmonic magnitudes should roughly decay: mag[k*f] < mag[(k-1)*f] */
    for (unsigned k = 2; k <= 5; k++)
        ok &= (mag[fb * k] < mag[fb * (k - 1)]);
    return ok;
}

/** Negative amplitude inverts the sawtooth wave. */
static int test_sawtooth_amplitude_negative(void)
{
    double sample_rate = 16000.0;
    double freq = 1000.0;
    unsigned period = (unsigned)(sample_rate / freq);
    unsigned N = period;
    double pos[N], neg[N];
    MD_sawtooth_wave(pos, N, 1.0, freq, sample_rate);
    MD_sawtooth_wave(neg, N, -1.0, freq, sample_rate);

    int ok = 1;
    for (unsigned i = 0; i < N; i++)
        ok &= approx_equal(neg[i], -pos[i], 1e-12);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_stft() and MD_stft_num_frames()
 * -----------------------------------------------------------------------*/

/** MD_stft_num_frames() formula and edge cases. */
static int test_stft_num_frames(void)
{
    int ok = 1;
    /* Standard case: 1024 samples, N=256, hop=128 -> (1024-256)/128+1 = 7 */
    ok &= (MD_stft_num_frames(1024, 256, 128) == 7);
    /* signal_len < N -> 0 frames */
    ok &= (MD_stft_num_frames(128, 256, 128) == 0);
    /* signal_len == 0 -> 0 frames */
    ok &= (MD_stft_num_frames(0, 256, 128) == 0);
    /* signal_len == N -> exactly 1 frame */
    ok &= (MD_stft_num_frames(256, 256, 128) == 1);
    /* hop == N, 4*N samples -> 4 frames */
    ok &= (MD_stft_num_frames(4 * 256, 256, 256) == 4);
    return ok;
}

/** All-zeros input should produce all-zero output. */
static int test_stft_silence(void)
{
    unsigned N = 256, hop = 128, signal_len = 1024;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    double *signal  = calloc(signal_len, sizeof(double));
    double *mag_out = malloc(num_frames * num_bins * sizeof(double));

    MD_stft(signal, signal_len, N, hop, mag_out);

    int ok = 1;
    for (unsigned i = 0; i < num_frames * num_bins; i++) {
        ok &= approx_equal(mag_out[i], 0.0, 1e-15);
    }

    free(mag_out);
    free(signal);
    return ok;
}

/** A bin-aligned pure tone should peak at the correct bin in every frame. */
static int test_stft_pure_tone(void)
{
    unsigned N           = 512;
    unsigned hop         = 256;
    double   sample_rate = 16000.0;
    double   freq        = 1000.0;   /* bin k0 = freq*N/sr = 32 */
    unsigned k0          = (unsigned)(freq * N / sample_rate);  /* = 32 */

    unsigned signal_len = 4 * N;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    double *signal  = malloc(signal_len * sizeof(double));
    double *mag_out = malloc(num_frames * num_bins * sizeof(double));

    /* Generate a pure tone that is bin-aligned */
    for (unsigned i = 0; i < signal_len; i++) {
        signal[i] = sin(2.0 * M_PI * freq * (double)i / sample_rate);
    }

    MD_stft(signal, signal_len, N, hop, mag_out);

    /* In every frame the peak bin (ignoring DC) should be k0 */
    int ok = 1;
    for (unsigned f = 0; f < num_frames; f++) {
        const double *row = mag_out + (size_t)f * num_bins;
        unsigned peak = 1;
        for (unsigned k = 2; k < num_bins; k++) {
            if (row[k] > row[peak]) peak = k;
        }
        ok &= (peak == k0);
    }

    free(mag_out);
    free(signal);
    return ok;
}

/** The number of output rows matches MD_stft_num_frames(), verified
 *  by checking the last row has non-trivial content. */
static int test_stft_frame_count(void)
{
    unsigned N = 256, hop = 128, signal_len = 2000;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    double *signal  = malloc(signal_len * sizeof(double));
    double *mag_out = calloc(num_frames * num_bins, sizeof(double));

    /* Use a sine so the output is non-trivial */
    for (unsigned i = 0; i < signal_len; i++) {
        signal[i] = sin(2.0 * M_PI * 10.0 * (double)i / (double)N);
    }

    MD_stft(signal, signal_len, N, hop, mag_out);

    /* The last row should have at least one non-zero magnitude */
    int ok = 0;
    const double *last_row = mag_out + (size_t)(num_frames - 1) * num_bins;
    for (unsigned k = 0; k < num_bins; k++) {
        if (last_row[k] > 1e-10) { ok = 1; break; }
    }

    free(mag_out);
    free(signal);
    return ok;
}

/** Non-overlapping frames (hop == N) from a 4N-sample signal. */
static int test_stft_hop_equal_N(void)
{
    unsigned N = 256, hop = 256, signal_len = 4 * N;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    /* Verify frame count */
    if (num_frames != 4) return 0;

    double *signal  = malloc(signal_len * sizeof(double));
    double *mag_out = malloc(num_frames * num_bins * sizeof(double));

    for (unsigned i = 0; i < signal_len; i++) {
        signal[i] = sin(2.0 * M_PI * 5.0 * (double)i / (double)N);
    }

    MD_stft(signal, signal_len, N, hop, mag_out);

    /* All magnitudes must be >= 0 */
    int ok = 1;
    for (unsigned i = 0; i < num_frames * num_bins; i++) {
        ok &= (mag_out[i] >= 0.0);
    }

    free(mag_out);
    free(signal);
    return ok;
}

/** All STFT magnitudes are >= 0 for a mixed-frequency input. */
static int test_stft_non_negative(void)
{
    unsigned N = 512, hop = 128, signal_len = 8192;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    double *signal  = malloc(signal_len * sizeof(double));
    double *mag_out = malloc(num_frames * num_bins * sizeof(double));

    for (unsigned i = 0; i < signal_len; i++) {
        signal[i] = 0.6 * sin(2.0 * M_PI * 200.0 * (double)i / 16000.0)
                  - 0.4 * cos(2.0 * M_PI * 800.0 * (double)i / 16000.0);
    }

    MD_stft(signal, signal_len, N, hop, mag_out);

    int ok = 1;
    for (unsigned i = 0; i < num_frames * num_bins; i++) {
        ok &= (mag_out[i] >= 0.0);
    }

    free(mag_out);
    free(signal);
    return ok;
}

/**
 * Parseval's theorem per frame:
 *   (mag[0]^2 + 2*sum(mag[1..N/2-1]^2) + mag[N/2]^2) / N
 *   == sum( (w[n]*x[n])^2 )
 *
 * This holds because MD_stft() returns unnormalised FFTW magnitudes,
 * consistent with MD_magnitude_spectrum().
 */
static int test_stft_parseval_per_frame(void)
{
    unsigned N = 512, hop = 256;
    unsigned signal_len = 4 * N;
    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    unsigned num_bins   = N / 2 + 1;

    double *signal  = malloc(signal_len * sizeof(double));
    double *window  = malloc(N * sizeof(double));
    double *mag_out = malloc(num_frames * num_bins * sizeof(double));

    MD_Gen_Hann_Win(window, N);

    for (unsigned i = 0; i < signal_len; i++) {
        signal[i] = 0.7 * sin(2.0 * M_PI * 50.0 * (double)i / (double)N)
                  + 0.3 * cos(2.0 * M_PI * 130.0 * (double)i / (double)N);
    }

    MD_stft(signal, signal_len, N, hop, mag_out);

    int ok = 1;
    for (unsigned f = 0; f < num_frames; f++) {
        const double *row = mag_out + (size_t)f * num_bins;
        const double *src = signal  + (size_t)f * hop;

        /* Time-domain energy after windowing */
        double td_energy = 0.0;
        for (unsigned n = 0; n < N; n++) {
            double s = src[n] * window[n];
            td_energy += s * s;
        }

        /* Frequency-domain energy via Parseval (unnormalised FFTW convention) */
        double fd_energy = row[0] * row[0];
        for (unsigned k = 1; k < N / 2; k++) {
            fd_energy += 2.0 * row[k] * row[k];
        }
        fd_energy += row[N / 2] * row[N / 2];
        fd_energy /= (double)N;

        ok &= approx_equal(td_energy, fd_energy, 1e-6);
    }

    free(mag_out);
    free(window);
    free(signal);
    return ok;
}

/** Calling with N=256 then N=512 should not crash and produce non-negative output. */
static int test_stft_plan_recache(void)
{
    int ok = 1;

    /* First call: N=256 with a sine wave (not impulse -- see plan notes) */
    {
        unsigned N = 256, hop = 64, signal_len = 4 * N;
        unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
        unsigned num_bins   = N / 2 + 1;

        double *signal  = malloc(signal_len * sizeof(double));
        double *mag_out = malloc(num_frames * num_bins * sizeof(double));

        for (unsigned i = 0; i < signal_len; i++) {
            signal[i] = sin(2.0 * M_PI * 10.0 * (double)i / (double)N);
        }

        MD_stft(signal, signal_len, N, hop, mag_out);

        for (unsigned i = 0; i < num_frames * num_bins; i++) {
            ok &= (mag_out[i] >= 0.0);
        }

        free(mag_out);
        free(signal);
    }

    /* Second call: N=512 (forces plan and window rebuild) */
    {
        unsigned N = 512, hop = 128, signal_len = 4 * N;
        unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
        unsigned num_bins   = N / 2 + 1;

        double *signal  = malloc(signal_len * sizeof(double));
        double *mag_out = malloc(num_frames * num_bins * sizeof(double));

        for (unsigned i = 0; i < signal_len; i++) {
            signal[i] = sin(2.0 * M_PI * 10.0 * (double)i / (double)N);
        }

        MD_stft(signal, signal_len, N, hop, mag_out);

        for (unsigned i = 0; i < num_frames * num_bins; i++) {
            ok &= (mag_out[i] >= 0.0);
        }

        free(mag_out);
        free(signal);
    }

    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for FIO_write_npy()
 * -----------------------------------------------------------------------*/

/** Write a 3x4 float matrix to .npy and verify the raw bytes. */
static int test_write_npy(void)
{
    const char *fname = "_test_output.npy";
    unlink(fname);

    /* Build a known 3x4 matrix */
    float row0[] = {1.0f, 2.0f, 3.0f, 4.0f};
    float row1[] = {5.0f, 6.0f, 7.0f, 8.0f};
    float row2[] = {9.0f, 10.0f, 11.0f, 12.0f};
    const float *rows[] = {row0, row1, row2};

    if (FIO_write_npy(fname, (const float **)rows, 3, 4) != 0) {
        unlink(fname);
        return 0;
    }

    /* Read back the file */
    FILE *f = fopen(fname, "rb");
    if (!f) { unlink(fname); return 0; }

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    unsigned char *buf = malloc((size_t)fsize);
    fread(buf, 1, (size_t)fsize, f);
    fclose(f);

    int ok = 1;

    /* Check magic bytes */
    ok &= (buf[0] == 0x93);
    ok &= (buf[1] == 'N' && buf[2] == 'U' && buf[3] == 'M' &&
            buf[4] == 'P' && buf[5] == 'Y');
    /* Check version 1.0 */
    ok &= (buf[6] == 1 && buf[7] == 0);

    /* Read header length (LE u16) */
    uint16_t hlen = (uint16_t)(buf[8] | (buf[9] << 8));

    /* Check total prefix is divisible by 64 */
    ok &= ((10 + hlen) % 64 == 0);

    /* Check header contains shape string */
    char *header = (char *)(buf + 10);
    ok &= (strstr(header, "'shape': (3, 4)") != nullptr);
    ok &= (strstr(header, "'descr': '<f4'") != nullptr);

    /* Check raw float data */
    size_t data_offset = 10 + hlen;
    float *data = (float *)(buf + data_offset);
    float expected[] = {1,2,3,4, 5,6,7,8, 9,10,11,12};
    for (int i = 0; i < 12; i++) {
        ok &= approx_equal((double)data[i], (double)expected[i], 1e-6);
    }

    /* Check file size */
    ok &= ((size_t)fsize == data_offset + 12 * sizeof(float));

    free(buf);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for FIO_write_safetensors()
 * -----------------------------------------------------------------------*/

/** Write known data to safetensors and verify the raw bytes. */
static int test_write_safetensors(void)
{
    const char *fname = "_test_output.safetensors";
    unlink(fname);

    float row0[] = {1.0f, 2.0f, 3.0f};
    float row1[] = {4.0f, 5.0f, 6.0f};
    const float *rows[] = {row0, row1};

    if (FIO_write_safetensors(fname, (const float **)rows, 2, 3) != 0) {
        unlink(fname);
        return 0;
    }

    FILE *f = fopen(fname, "rb");
    if (!f) { unlink(fname); return 0; }

    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);

    unsigned char *buf = malloc((size_t)fsize);
    fread(buf, 1, (size_t)fsize, f);
    fclose(f);

    int ok = 1;

    /* Read 8-byte LE u64 header size */
    uint64_t hsize = 0;
    for (int i = 0; i < 8; i++) {
        hsize |= ((uint64_t)buf[i]) << (i * 8);
    }

    /* JSON header should be reasonable size */
    ok &= (hsize > 10 && hsize < 500);

    /* Check JSON contains expected fields */
    char *json = malloc((size_t)hsize + 1);
    memcpy(json, buf + 8, (size_t)hsize);
    json[hsize] = '\0';

    ok &= (strstr(json, "\"dtype\":\"F32\"") != nullptr);
    ok &= (strstr(json, "\"shape\":[2,3]") != nullptr);
    ok &= (strstr(json, "\"data_offsets\":[0,24]") != nullptr);  /* 2*3*4 = 24 */

    free(json);

    /* Check raw float data */
    size_t data_offset = 8 + (size_t)hsize;
    float *data = (float *)(buf + data_offset);
    float expected[] = {1,2,3, 4,5,6};
    for (int i = 0; i < 6; i++) {
        ok &= approx_equal((double)data[i], (double)expected[i], 1e-6);
    }

    /* Check file size */
    ok &= ((size_t)fsize == data_offset + 6 * sizeof(float));

    free(buf);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for FIO_write_wav()
 * -----------------------------------------------------------------------*/

/** Write a sine wave to WAV, read it back, verify round-trip. */
static int test_write_wav(void)
{
    const char *fname = "_test_output.wav";
    unlink(fname);

    unsigned samprate = 16000;
    size_t datalen = 1600;  /* 0.1 seconds */
    float *data = malloc(datalen * sizeof(float));

    /* Generate a 440 Hz sine wave */
    for (size_t i = 0; i < datalen; i++) {
        data[i] = (float)sin(2.0 * M_PI * 440.0 * (double)i / (double)samprate);
    }

    if (FIO_write_wav(fname, data, datalen, samprate) != 0) {
        free(data);
        unlink(fname);
        return 0;
    }

    /* Read it back */
    float *readback = nullptr;
    size_t readlen = 0;
    unsigned readrate = 0;
    if (FIO_read_audio(fname, &readback, &readlen, &readrate, 0) != 0) {
        free(data);
        unlink(fname);
        return 0;
    }

    int ok = 1;

    /* Check sample count and rate */
    ok &= (readlen == datalen);
    ok &= (readrate == samprate);

    /* Check values round-trip within tolerance */
    if (readlen == datalen) {
        for (size_t i = 0; i < datalen; i++) {
            ok &= approx_equal((double)data[i], (double)readback[i], 1e-5);
        }
    }

    free(readback);
    free(data);
    unlink(fname);
    return ok;
}

/* -----------------------------------------------------------------------
 * Main: run all tests
 * -----------------------------------------------------------------------*/

int main(void)
{
    printf("=== miniDSP Test Suite ===\n\n");

    printf("--- MD_dot ---\n");
    RUN_TEST(test_dot_orthogonal);
    RUN_TEST(test_dot_self);
    RUN_TEST(test_dot_known);
    RUN_TEST(test_dot_single);
    RUN_TEST(test_dot_zeros);

    printf("\n--- MD_energy ---\n");
    RUN_TEST(test_energy_zero);
    RUN_TEST(test_energy_constant);
    RUN_TEST(test_energy_known);
    RUN_TEST(test_energy_single);

    printf("\n--- MD_power ---\n");
    RUN_TEST(test_power_known);
    RUN_TEST(test_power_sine);

    printf("\n--- MD_power_db ---\n");
    RUN_TEST(test_power_db_known);
    RUN_TEST(test_power_db_quiet);
    RUN_TEST(test_power_db_floor);

    printf("\n--- MD_scale / MD_scale_vec ---\n");
    RUN_TEST(test_scale_midpoint);
    RUN_TEST(test_scale_endpoints);
    RUN_TEST(test_scale_vec);

    printf("\n--- MD_fit_within_range ---\n");
    RUN_TEST(test_fit_in_range_no_change);
    RUN_TEST(test_fit_in_range_rescale);

    printf("\n--- MD_adjust_dblevel ---\n");
    RUN_TEST(test_adjust_dblevel);

    printf("\n--- MD_entropy ---\n");
    RUN_TEST(test_entropy_uniform);
    RUN_TEST(test_entropy_spike);
    RUN_TEST(test_entropy_zeros);
    RUN_TEST(test_entropy_single);
    RUN_TEST(test_entropy_no_clip);
    RUN_TEST(test_entropy_clip);

    printf("\n--- MD_rms ---\n");
    RUN_TEST(test_rms_sine);
    RUN_TEST(test_rms_dc);
    RUN_TEST(test_rms_silence);
    RUN_TEST(test_rms_matches_sqrt_power);

    printf("\n--- MD_zero_crossing_rate ---\n");
    RUN_TEST(test_zcr_sine);
    RUN_TEST(test_zcr_constant);
    RUN_TEST(test_zcr_alternating);
    RUN_TEST(test_zcr_noise_higher_than_sine);
    RUN_TEST(test_zcr_zeros_handling);

    printf("\n--- MD_autocorrelation ---\n");
    RUN_TEST(test_acf_lag0);
    RUN_TEST(test_acf_sine_period);
    RUN_TEST(test_acf_noise_decay);
    RUN_TEST(test_acf_bounded);
    RUN_TEST(test_acf_silence);

    printf("\n--- MD_peak_detect ---\n");
    RUN_TEST(test_peaks_known);
    RUN_TEST(test_peaks_threshold);
    RUN_TEST(test_peaks_min_distance);
    RUN_TEST(test_peaks_none);
    RUN_TEST(test_peaks_flat_signal);
    RUN_TEST(test_peaks_single_sample);

    printf("\n--- Pitch detection ---\n");
    RUN_TEST(test_f0_acf_clean_sine);
    RUN_TEST(test_f0_acf_harmonic_signal);
    RUN_TEST(test_f0_acf_noisy_sine);
    RUN_TEST(test_f0_acf_silence);
    RUN_TEST(test_f0_acf_out_of_range);
    RUN_TEST(test_f0_acf_lag_edge_mapping);
    RUN_TEST(test_f0_fft_clean_sine);
    RUN_TEST(test_f0_fft_two_tone_dominance);
    RUN_TEST(test_f0_fft_silence);
    RUN_TEST(test_f0_fft_out_of_range);
    RUN_TEST(test_f0_fft_bin_edge_mapping);
    RUN_TEST(test_f0_fft_plan_recache);

    printf("\n--- MD_mix ---\n");
    RUN_TEST(test_mix_equal_weights);
    RUN_TEST(test_mix_passthrough);
    RUN_TEST(test_mix_energy);
    RUN_TEST(test_mix_inplace);

    printf("\n--- Simple effects ---\n");
    RUN_TEST(test_delay_echo_dry_passthrough);
    RUN_TEST(test_delay_echo_impulse_decay);
    RUN_TEST(test_delay_echo_inplace);
    RUN_TEST(test_tremolo_depth_zero_passthrough);
    RUN_TEST(test_tremolo_gain_bounds);
    RUN_TEST(test_comb_reverb_dry_passthrough);
    RUN_TEST(test_comb_reverb_impulse_decay);
    RUN_TEST(test_comb_reverb_inplace);

    printf("\n--- FIR filters / convolution ---\n");
    RUN_TEST(test_convolution_num_samples);
    RUN_TEST(test_convolution_time_known_small);
    RUN_TEST(test_convolution_time_impulse_identity);
    RUN_TEST(test_fir_filter_known_taps);
    RUN_TEST(test_moving_average_step_response);
    RUN_TEST(test_moving_average_matches_boxcar_fir);
    RUN_TEST(test_convolution_fft_ola_matches_time);
    RUN_TEST(test_convolution_fft_ola_different_lengths);

    printf("\n--- Window generation ---\n");
    RUN_TEST(test_hann_endpoints);
    RUN_TEST(test_hann_peak);
    RUN_TEST(test_hann_symmetry);
    RUN_TEST(test_hann_range);
    RUN_TEST(test_hamming_endpoints);
    RUN_TEST(test_hamming_peak);
    RUN_TEST(test_hamming_symmetry);
    RUN_TEST(test_hamming_range);
    RUN_TEST(test_blackman_endpoints);
    RUN_TEST(test_blackman_symmetry);
    RUN_TEST(test_blackman_range);
    RUN_TEST(test_rect_all_ones);
    RUN_TEST(test_window_single_sample);

    printf("\n--- MD_magnitude_spectrum ---\n");
    RUN_TEST(test_mag_spectrum_single_sine);
    RUN_TEST(test_mag_spectrum_two_sines);
    RUN_TEST(test_mag_spectrum_dc);
    RUN_TEST(test_mag_spectrum_zeros);
    RUN_TEST(test_mag_spectrum_impulse);
    RUN_TEST(test_mag_spectrum_parseval);
    RUN_TEST(test_mag_spectrum_different_lengths);
    RUN_TEST(test_mag_spectrum_non_negative);

    printf("\n--- MD_power_spectral_density ---\n");
    RUN_TEST(test_psd_single_sine);
    RUN_TEST(test_psd_two_sines);
    RUN_TEST(test_psd_dc);
    RUN_TEST(test_psd_zeros);
    RUN_TEST(test_psd_impulse);
    RUN_TEST(test_psd_parseval);
    RUN_TEST(test_psd_different_lengths);
    RUN_TEST(test_psd_non_negative);

    printf("\n--- MD_phase_spectrum ---\n");
    RUN_TEST(test_phase_spectrum_cosine);
    RUN_TEST(test_phase_spectrum_sine);
    RUN_TEST(test_phase_spectrum_impulse_at_zero);
    RUN_TEST(test_phase_spectrum_impulse_at_one);
    RUN_TEST(test_phase_spectrum_zeros);
    RUN_TEST(test_phase_spectrum_different_lengths);

    printf("\n--- GCC-PHAT delay estimation ---\n");
    RUN_TEST(test_gcc_phat_positive_delay);
    RUN_TEST(test_gcc_phat_negative_delay);
    RUN_TEST(test_gcc_phat_zero_delay);
    RUN_TEST(test_gcc_simp_delay);
    RUN_TEST(test_gcc_returns_entropy);
    RUN_TEST(test_gcc_multiple_delays);
    RUN_TEST(test_gcc_lagvals_structure);
    RUN_TEST(test_gcc_different_lengths);

    printf("\n--- Biquad filters ---\n");
    RUN_TEST(test_biquad_invalid_type);
    RUN_TEST(test_biquad_lpf);
    RUN_TEST(test_biquad_hpf);
    RUN_TEST(test_biquad_bpf);
    RUN_TEST(test_biquad_notch);
    RUN_TEST(test_biquad_peq_boost);
    RUN_TEST(test_biquad_low_shelf);
    RUN_TEST(test_biquad_high_shelf);
    RUN_TEST(test_biquad_hpf_dc_rejection);

    printf("\n--- MD_sine_wave ---\n");
    RUN_TEST(test_sine_starts_at_zero);
    RUN_TEST(test_sine_quarter_period);
    RUN_TEST(test_sine_full_period);
    RUN_TEST(test_sine_spectrum_peak);
    RUN_TEST(test_sine_amplitude_negative);

    printf("\n--- MD_white_noise ---\n");
    RUN_TEST(test_white_noise_mean_near_zero);
    RUN_TEST(test_white_noise_stddev_matches_amplitude);
    RUN_TEST(test_white_noise_reproducible);
    RUN_TEST(test_white_noise_different_seeds);
    RUN_TEST(test_white_noise_flat_spectrum);
    RUN_TEST(test_white_noise_odd_length);

    printf("\n--- MD_impulse ---\n");
    RUN_TEST(test_impulse_at_zero);
    RUN_TEST(test_impulse_at_middle);
    RUN_TEST(test_impulse_at_last);
    RUN_TEST(test_impulse_amplitude);
    RUN_TEST(test_impulse_amplitude_negative);
    RUN_TEST(test_impulse_single_sample);
    RUN_TEST(test_impulse_flat_spectrum);

    printf("\n--- MD_chirp_linear ---\n");
    RUN_TEST(test_chirp_linear_starts_at_zero);
    RUN_TEST(test_chirp_linear_spectrum_spread);
    RUN_TEST(test_chirp_linear_constant_freq);
    RUN_TEST(test_chirp_linear_amplitude_bound);

    printf("\n--- MD_chirp_log ---\n");
    RUN_TEST(test_chirp_log_starts_at_zero);
    RUN_TEST(test_chirp_log_spectrum_spread);
    RUN_TEST(test_chirp_log_amplitude_bound);

    printf("\n--- MD_square_wave ---\n");
    RUN_TEST(test_square_high_low);
    RUN_TEST(test_square_spectrum_harmonics);
    RUN_TEST(test_square_amplitude_negative);

    printf("\n--- MD_sawtooth_wave ---\n");
    RUN_TEST(test_sawtooth_midpoint);
    RUN_TEST(test_sawtooth_spectrum_harmonics);
    RUN_TEST(test_sawtooth_amplitude_negative);

    printf("\n--- MD_stft ---\n");
    RUN_TEST(test_stft_num_frames);
    RUN_TEST(test_stft_silence);
    RUN_TEST(test_stft_pure_tone);
    RUN_TEST(test_stft_frame_count);
    RUN_TEST(test_stft_hop_equal_N);
    RUN_TEST(test_stft_non_negative);
    RUN_TEST(test_stft_parseval_per_frame);
    RUN_TEST(test_stft_plan_recache);

    printf("\n--- File I/O writers ---\n");
    RUN_TEST(test_write_npy);
    RUN_TEST(test_write_safetensors);
    RUN_TEST(test_write_wav);

    /* Clean up FFTW resources */
    MD_shutdown();

    /* Print summary */
    printf("\n=== Results: %d/%d passed", tests_passed, tests_run);
    if (tests_failed > 0) {
        printf(", %d FAILED", tests_failed);
    }
    printf(" ===\n");

    return (tests_failed > 0) ? 1 : 0;
}
