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
#include "../fileio.h"

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
 * Tests for MD_Gen_Hann_Win()
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

    printf("\n--- MD_Gen_Hann_Win ---\n");
    RUN_TEST(test_hann_endpoints);
    RUN_TEST(test_hann_peak);
    RUN_TEST(test_hann_symmetry);
    RUN_TEST(test_hann_range);

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
