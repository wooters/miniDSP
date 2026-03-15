/**
 * @file test_core.c
 * @brief Tests for core DSP functions in minidsp_core.c.
 *
 * Covers: MD_dot, MD_energy, MD_power, MD_power_db, MD_scale, MD_scale_vec,
 * MD_fit_within_range, MD_adjust_dblevel, MD_entropy, MD_rms,
 * MD_zero_crossing_rate, MD_autocorrelation, MD_peak_detect,
 * MD_f0_autocorrelation, MD_f0_fft, MD_mix.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

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

/** Power of an all-zero signal should be zero. */
static int test_power_zero(void)
{
    double a[] = {0.0, 0.0, 0.0, 0.0};
    return approx_equal(MD_power(a, 4), 0.0, 1e-15);
}

/** Power of a single sample equals x^2. */
static int test_power_single(void)
{
    double a[] = {5.0};
    return approx_equal(MD_power(a, 1), 25.0, 1e-15);
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

/** Scaling with identical in/out ranges returns input unchanged. */
static int test_scale_identity_range(void)
{
    int ok = 1;
    ok &= approx_equal(MD_scale(3.7, 0.0, 10.0, 0.0, 10.0), 3.7, 1e-15);
    ok &= approx_equal(MD_scale(-5.0, -10.0, 10.0, -10.0, 10.0), -5.0, 1e-15);
    return ok;
}

/** Inverted output range reverses direction. */
static int test_scale_inverted_range(void)
{
    int ok = 1;
    ok &= approx_equal(MD_scale(0.0, 0.0, 10.0, 100.0, 0.0), 100.0, 1e-15);
    ok &= approx_equal(MD_scale(10.0, 0.0, 10.0, 100.0, 0.0), 0.0, 1e-15);
    ok &= approx_equal(MD_scale(5.0, 0.0, 10.0, 100.0, 0.0), 50.0, 1e-15);
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

/** Constant signal (min==max): already fits, copied unchanged. */
static int test_fit_in_range_constant(void)
{
    double in[]  = {7.0, 7.0, 7.0};
    double out[3];
    MD_fit_within_range(in, out, 3, 0.0, 10.0);
    int ok = 1;
    ok &= approx_equal(out[0], 7.0, 1e-15);
    ok &= approx_equal(out[1], 7.0, 1e-15);
    ok &= approx_equal(out[2], 7.0, 1e-15);
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

/** Silent input should produce silent output (not NaN/Inf). */
static int test_adjust_dblevel_silence(void)
{
    unsigned N = 100;
    double in[100], out[100];
    memset(in, 0, sizeof(in));

    MD_adjust_dblevel(in, out, N, -10.0);

    for (unsigned i = 0; i < N; i++) {
        if (out[i] != 0.0) return 0;
    }
    return 1;
}

/** Extremely quiet input should produce finite output (no NaN/Inf). */
static int test_adjust_dblevel_near_zero(void)
{
    unsigned N = 100;
    double in[100], out[100];
    for (unsigned i = 0; i < N; i++) in[i] = 1e-300;

    MD_adjust_dblevel(in, out, N, -10.0);

    for (unsigned i = 0; i < N; i++) {
        if (isnan(out[i]) || isinf(out[i])) return 0;
    }
    return 1;
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

/** Both weights zero should produce silence. */
static int test_mix_zero_weights(void)
{
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {4.0, 5.0, 6.0};
    double out[3];
    MD_mix(a, b, out, 3, 0.0, 0.0);
    int ok = 1;
    ok &= approx_equal(out[0], 0.0, 1e-15);
    ok &= approx_equal(out[1], 0.0, 1e-15);
    ok &= approx_equal(out[2], 0.0, 1e-15);
    return ok;
}

/** Weight 0.0/1.0 should pass through the second signal. */
static int test_mix_passthrough_b(void)
{
    double a[] = {1.0, 2.0, 3.0};
    double b[] = {10.0, 20.0, 30.0};
    double out[3];
    MD_mix(a, b, out, 3, 0.0, 1.0);
    int ok = 1;
    ok &= approx_equal(out[0], 10.0, 1e-15);
    ok &= approx_equal(out[1], 20.0, 1e-15);
    ok &= approx_equal(out[2], 30.0, 1e-15);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public runner
 * -----------------------------------------------------------------------*/

void run_core_tests(void)
{
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
    RUN_TEST(test_power_zero);
    RUN_TEST(test_power_single);

    printf("\n--- MD_power_db ---\n");
    RUN_TEST(test_power_db_known);
    RUN_TEST(test_power_db_quiet);
    RUN_TEST(test_power_db_floor);

    printf("\n--- MD_scale / MD_scale_vec ---\n");
    RUN_TEST(test_scale_midpoint);
    RUN_TEST(test_scale_endpoints);
    RUN_TEST(test_scale_identity_range);
    RUN_TEST(test_scale_inverted_range);
    RUN_TEST(test_scale_vec);

    printf("\n--- MD_fit_within_range ---\n");
    RUN_TEST(test_fit_in_range_no_change);
    RUN_TEST(test_fit_in_range_rescale);
    RUN_TEST(test_fit_in_range_constant);

    printf("\n--- MD_adjust_dblevel ---\n");
    RUN_TEST(test_adjust_dblevel);
    RUN_TEST(test_adjust_dblevel_silence);
    RUN_TEST(test_adjust_dblevel_near_zero);

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
    RUN_TEST(test_mix_zero_weights);
    RUN_TEST(test_mix_passthrough_b);
}
