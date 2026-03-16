/**
 * @file test_spectrum.c
 * @brief Tests for spectrum analysis: magnitude, PSD, phase, STFT, mel, MFCC.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

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

/** Phase of a constant (DC-only) signal: DC bin should have zero phase. */
static int test_phase_spectrum_dc_signal(void)
{
    unsigned N = 256;
    double sig[256];
    for (unsigned i = 0; i < N; i++) sig[i] = 3.0;

    unsigned num_bins = N / 2 + 1;
    double phase[num_bins];
    MD_phase_spectrum(sig, N, phase);

    /* DC bin should be 0 phase (all energy is real, positive) */
    int ok = approx_equal(phase[0], 0.0, 1e-10);
    /* All non-DC bins should also have ~0 phase (they're ~0 magnitude) */
    for (unsigned k = 1; k < num_bins; k++) {
        ok &= approx_equal(phase[k], 0.0, 1e-10);
    }
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

/** Single-frame STFT (signal_len == N) should match windowed magnitude spectrum. */
static int test_stft_single_frame(void)
{
    unsigned N = 256;
    unsigned num_bins = N / 2 + 1;
    double signal[256];

    /* Create a simple signal */
    MD_sine_wave(signal, N, 1.0, 1000.0, 8000.0);

    /* STFT with signal_len == N -> exactly 1 frame */
    unsigned num_frames = MD_stft_num_frames(N, N, N);
    if (num_frames != 1) return 0;

    double stft_out[num_bins];
    MD_stft(signal, N, N, N, stft_out);

    /* Compute magnitude spectrum of the Hann-windowed signal for comparison */
    double windowed[256];
    double window[256];
    MD_Gen_Hann_Win(window, N);
    for (unsigned i = 0; i < N; i++) windowed[i] = signal[i] * window[i];

    double mag[num_bins];
    MD_magnitude_spectrum(windowed, N, mag);

    /* They should match */
    int ok = 1;
    for (unsigned k = 0; k < num_bins; k++) {
        ok &= approx_equal(stft_out[k], mag[k], 1e-10);
    }
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for mel filterbanks and MFCCs
 * -----------------------------------------------------------------------*/

/** Mel filterbank weights should be finite, non-negative, bounded, and triangular. */
static int test_mel_filterbank_bounds_and_triangles(void)
{
    const unsigned N = 512;
    const unsigned num_mels = 20;
    const double sample_rate = 16000.0;
    const unsigned num_bins = N / 2 + 1;

    double *fb = malloc((size_t)num_mels * num_bins * sizeof(double));
    if (!fb) return 0;

    MD_mel_filterbank(N, sample_rate, num_mels, 80.0, 7600.0, fb);

    int ok = 1;
    for (unsigned m = 0; m < num_mels; m++) {
        const double *row = fb + (size_t)m * num_bins;

        unsigned peak = 0;
        for (unsigned k = 0; k < num_bins; k++) {
            double w = row[k];
            ok &= isfinite(w);
            ok &= (w >= -1e-12);
            ok &= (w <= 1.0 + 1e-12);
            if (w > row[peak]) peak = k;
        }

        if (row[peak] <= 1e-12) continue;

        for (unsigned k = 1; k <= peak; k++) {
            ok &= (row[k] + 1e-10 >= row[k - 1]);
        }
        for (unsigned k = peak + 1; k < num_bins; k++) {
            ok &= (row[k] <= row[k - 1] + 1e-10);
        }
    }

    free(fb);
    return ok;
}

/** A high-frequency tone should peak at a higher mel band than a low-frequency tone. */
static int test_mel_energies_tone_band_shift(void)
{
    const unsigned N = 1024;
    const unsigned num_mels = 26;
    const double sample_rate = 16000.0;

    double low[N], high[N];
    double mel_low[num_mels], mel_high[num_mels];

    MD_sine_wave(low, N, 1.0, 200.0, sample_rate);
    MD_sine_wave(high, N, 1.0, 3200.0, sample_rate);

    MD_mel_energies(low, N, sample_rate, num_mels, 80.0, 7600.0, mel_low);
    MD_mel_energies(high, N, sample_rate, num_mels, 80.0, 7600.0, mel_high);

    unsigned peak_low = 0, peak_high = 0;
    for (unsigned m = 1; m < num_mels; m++) {
        if (mel_low[m] > mel_low[peak_low]) peak_low = m;
        if (mel_high[m] > mel_high[peak_high]) peak_high = m;
    }

    return peak_high > peak_low;
}

/** Silence should produce all-zero mel energies. */
static int test_mel_energies_silence(void)
{
    const unsigned N = 1024;
    const unsigned num_mels = 26;
    double sig[N];
    double mel[num_mels];

    memset(sig, 0, sizeof(sig));
    MD_mel_energies(sig, N, 16000.0, num_mels, 80.0, 7600.0, mel);

    int ok = 1;
    for (unsigned m = 0; m < num_mels; m++) {
        ok &= approx_equal(mel[m], 0.0, 1e-15);
    }
    return ok;
}

/** Silence MFCCs should be finite; higher coefficients should be near zero. */
static int test_mfcc_silence_finite(void)
{
    const unsigned N = 1024;
    const unsigned num_mels = 26;
    const unsigned num_coeffs = 13;
    double sig[N];
    double mfcc[num_coeffs];

    memset(sig, 0, sizeof(sig));
    MD_mfcc(sig, N, 16000.0, num_mels, num_coeffs, 80.0, 7600.0, mfcc);

    int ok = 1;
    for (unsigned c = 0; c < num_coeffs; c++) {
        ok &= isfinite(mfcc[c]);
    }
    for (unsigned c = 1; c < num_coeffs; c++) {
        ok &= (fabs(mfcc[c]) < 1e-10);
    }
    return ok;
}

/** Scaling a frame should mostly change C0; higher MFCCs should remain similar. */
static int test_mfcc_amplitude_scaling_invariance(void)
{
    const unsigned N = 1024;
    const unsigned num_mels = 26;
    const unsigned num_coeffs = 13;
    const double scale = 0.4;

    double noise[N], scaled[N];
    double mfcc_ref[num_coeffs], mfcc_scaled[num_coeffs];

    MD_white_noise(noise, N, 0.4, 42);
    for (unsigned i = 0; i < N; i++) {
        scaled[i] = scale * noise[i];
    }

    MD_mfcc(noise, N, 16000.0, num_mels, num_coeffs, 80.0, 7600.0, mfcc_ref);
    MD_mfcc(scaled, N, 16000.0, num_mels, num_coeffs, 80.0, 7600.0, mfcc_scaled);

    int ok = 1;
    for (unsigned c = 1; c < num_coeffs; c++) {
        ok &= approx_equal(mfcc_ref[c], mfcc_scaled[c], 1e-3);
    }

    /* Power-domain scaling by scale^2 adds a constant in log-mel space,
     * which maps only to C0 under orthonormal DCT-II. */
    double expected_delta_c0 = sqrt((double)num_mels) * log(scale * scale);
    double delta_c0 = mfcc_scaled[0] - mfcc_ref[0];
    ok &= approx_equal(delta_c0, expected_delta_c0, 1e-3);

    return ok;
}

/** Consecutive calls with different frame sizes should rebuild caches safely. */
static int test_mel_mfcc_plan_recache(void)
{
    int ok = 1;

    {
        const unsigned N = 1024;
        const unsigned num_mels = 26;
        double sig[N], mel[num_mels], mfcc[13];
        MD_sine_wave(sig, N, 1.0, 440.0, 16000.0);
        MD_mel_energies(sig, N, 16000.0, num_mels, 80.0, 7600.0, mel);
        MD_mfcc(sig, N, 16000.0, num_mels, 13, 80.0, 7600.0, mfcc);
        for (unsigned i = 0; i < num_mels; i++) ok &= isfinite(mel[i]);
        for (unsigned i = 0; i < 13; i++) ok &= isfinite(mfcc[i]);
    }

    {
        const unsigned N = 512;
        const unsigned num_mels = 20;
        double sig[N], mel[num_mels], mfcc[12];
        MD_sine_wave(sig, N, 1.0, 880.0, 8000.0);
        MD_mel_energies(sig, N, 8000.0, num_mels, 50.0, 3900.0, mel);
        MD_mfcc(sig, N, 8000.0, num_mels, 12, 50.0, 3900.0, mfcc);
        for (unsigned i = 0; i < num_mels; i++) ok &= isfinite(mel[i]);
        for (unsigned i = 0; i < 12; i++) ok &= isfinite(mfcc[i]);
    }

    return ok;
}

/** Out-of-range frequency bounds should clamp to [0, Nyquist]. */
static int test_mel_frequency_clamp_behavior(void)
{
    const unsigned N = 1024;
    const unsigned num_mels = 26;
    const double sample_rate = 16000.0;
    double sig[N];
    double mel_clamped[num_mels], mel_ref[num_mels], mel_empty[num_mels];
    double mfcc_empty[13];

    MD_sine_wave(sig, N, 1.0, 440.0, sample_rate);

    MD_mel_energies(sig, N, sample_rate, num_mels, -500.0, 50000.0, mel_clamped);
    MD_mel_energies(sig, N, sample_rate, num_mels, 0.0, sample_rate / 2.0, mel_ref);

    int ok = 1;
    for (unsigned m = 0; m < num_mels; m++) {
        ok &= approx_equal(mel_clamped[m], mel_ref[m], 1e-10);
    }

    /* Empty clamped band (both edges at Nyquist) should yield all-zero mel energies. */
    MD_mel_energies(sig, N, sample_rate, num_mels, 9000.0, 10000.0, mel_empty);
    for (unsigned m = 0; m < num_mels; m++) {
        ok &= approx_equal(mel_empty[m], 0.0, 1e-15);
    }

    /* MFCCs should still be finite due log floor handling. */
    MD_mfcc(sig, N, sample_rate, num_mels, 13, 9000.0, 10000.0, mfcc_empty);
    for (unsigned c = 0; c < 13; c++) {
        ok &= isfinite(mfcc_empty[c]);
    }

    return ok;
}

/** Extremely narrow mel ranges should produce degenerate all-zero rows safely. */
static int test_mel_filterbank_degenerate_geometry(void)
{
    const unsigned N = 16;
    const unsigned num_mels = 40;
    const unsigned num_bins = N / 2 + 1;
    double fb[(size_t)num_mels * num_bins];

    MD_mel_filterbank(N, 16000.0, num_mels, 1000.0, 1010.0, fb);

    int ok = 1;
    unsigned zero_rows = 0;
    for (unsigned m = 0; m < num_mels; m++) {
        const double *row = fb + (size_t)m * num_bins;
        double row_sum = 0.0;
        for (unsigned k = 0; k < num_bins; k++) {
            ok &= isfinite(row[k]);
            ok &= (row[k] >= -1e-12);
            row_sum += row[k];
        }
        if (fabs(row_sum) < 1e-12) zero_rows++;
    }

    ok &= (zero_rows > 0);
    return ok;
}

/** Deterministic MFCC regression for a fixed frame and parameter set.
 *
 * Golden values were generated with an independent Python direct-DFT script using:
 * - HTK mel mapping
 * - Hann window
 * - one-sided PSD (|X|^2 / N)
 * - ln(max(E, 1e-12))
 * - orthonormal DCT-II with C0 included */
static int test_mfcc_golden_reference(void)
{
    const unsigned N = 512;
    const unsigned num_mels = 26;
    const unsigned num_coeffs = 13;
    const double fs = 16000.0;
    double sig[N];
    double mfcc[num_coeffs];

    for (unsigned n = 0; n < N; n++) {
        double t = (double)n / fs;
        sig[n] = 0.7 * sin(2.0 * M_PI * 440.0 * t)
               + 0.2 * cos(2.0 * M_PI * 1000.0 * t)
               + 0.1 * sin(2.0 * M_PI * 3000.0 * t);
    }

    MD_mfcc(sig, N, fs, num_mels, num_coeffs, 80.0, 7600.0, mfcc);

    const double expected[13] = {
        -75.381495813245806, 36.882955249179645, -4.552916285440863,
        3.775663546970039, -18.417919467564612, -10.745741736169254,
        12.729163817453310, -7.496995440686327, -17.500169058663040,
        -2.934673446412961, -8.198452514267208, -3.352676395931926,
        15.283540889806705
    };

    int ok = 1;
    for (unsigned c = 0; c < num_coeffs; c++) {
        ok &= approx_equal(mfcc[c], expected[c], 1e-6);
    }
    return ok;
}

/* -----------------------------------------------------------------------
 * MD_lowpass_brickwall tests
 * -----------------------------------------------------------------------*/

/** A tone below the cutoff should be preserved (energy unchanged). */
static int test_lowpass_brickwall_preserves_low_freq(void)
{
    unsigned N = 48000;
    double sr = 48000.0;
    double freq = 1000.0;
    double *sig = malloc(N * sizeof(double));
    MD_sine_wave(sig, N, 1.0, freq, sr);

    double energy_before = MD_energy(sig, N);
    MD_lowpass_brickwall(sig, N, 8000.0, sr);
    double energy_after = MD_energy(sig, N);

    free(sig);
    /* Within 0.1% */
    return approx_equal(energy_after, energy_before,
                        energy_before * 0.001);
}

/** A tone above the cutoff should be eliminated (energy ~ 0). */
static int test_lowpass_brickwall_removes_high_freq(void)
{
    unsigned N = 48000;
    double sr = 48000.0;
    double freq = 20000.0;
    double *sig = malloc(N * sizeof(double));
    MD_sine_wave(sig, N, 1.0, freq, sr);

    MD_lowpass_brickwall(sig, N, 8000.0, sr);
    double energy_after = MD_energy(sig, N);

    free(sig);
    return (energy_after < 1e-18);
}

/** Mixed signal: only the low-frequency component should survive. */
static int test_lowpass_brickwall_mixed_signal(void)
{
    unsigned N = 48000;
    double sr = 48000.0;
    double *lo = malloc(N * sizeof(double));
    double *sig = malloc(N * sizeof(double));
    MD_sine_wave(lo, N, 1.0, 1000.0, sr);

    /* Build mixed signal: 1000 Hz + 20000 Hz */
    for (unsigned i = 0; i < N; i++)
        sig[i] = lo[i] + sin(2.0 * M_PI * 20000.0 * (double)i / sr);

    double energy_lo = MD_energy(lo, N);
    MD_lowpass_brickwall(sig, N, 8000.0, sr);
    double energy_after = MD_energy(sig, N);

    /* After filtering, energy should match the low-freq component */
    int ok = approx_equal(energy_after, energy_lo, energy_lo * 0.001);

    /* Verify the high-freq residual is gone by checking difference */
    double diff_energy = 0.0;
    for (unsigned i = 0; i < N; i++) {
        double d = sig[i] - lo[i];
        diff_energy += d * d;
    }
    ok &= (diff_energy < 1e-14);

    free(sig);
    free(lo);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_spectrum_tests(void)
{
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
    RUN_TEST(test_phase_spectrum_dc_signal);

    printf("\n--- MD_stft ---\n");
    RUN_TEST(test_stft_num_frames);
    RUN_TEST(test_stft_silence);
    RUN_TEST(test_stft_pure_tone);
    RUN_TEST(test_stft_frame_count);
    RUN_TEST(test_stft_hop_equal_N);
    RUN_TEST(test_stft_non_negative);
    RUN_TEST(test_stft_parseval_per_frame);
    RUN_TEST(test_stft_plan_recache);
    RUN_TEST(test_stft_single_frame);

    printf("\n--- Mel / MFCC ---\n");
    RUN_TEST(test_mel_filterbank_bounds_and_triangles);
    RUN_TEST(test_mel_energies_tone_band_shift);
    RUN_TEST(test_mel_energies_silence);
    RUN_TEST(test_mfcc_silence_finite);
    RUN_TEST(test_mfcc_amplitude_scaling_invariance);
    RUN_TEST(test_mel_mfcc_plan_recache);
    RUN_TEST(test_mel_frequency_clamp_behavior);
    RUN_TEST(test_mel_filterbank_degenerate_geometry);
    RUN_TEST(test_mfcc_golden_reference);

    printf("\n--- MD_lowpass_brickwall ---\n");
    RUN_TEST(test_lowpass_brickwall_preserves_low_freq);
    RUN_TEST(test_lowpass_brickwall_removes_high_freq);
    RUN_TEST(test_lowpass_brickwall_mixed_signal);
}
