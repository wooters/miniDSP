/**
 * @file test_generators.c
 * @brief Tests for signal generator functions.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "minidsp.h"
#include "test_helpers.h"

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

/** Zero amplitude produces all-zero output. */
static int test_sine_zero_amplitude(void)
{
    double out[64];
    MD_sine_wave(out, 64, 0.0, 440.0, 44100.0);
    for (unsigned i = 0; i < 64; i++) {
        if (fabs(out[i]) > 1e-15) return 0;
    }
    return 1;
}

/** Zero frequency: sin(0) = 0 for all samples. */
static int test_sine_zero_frequency(void)
{
    double out[64];
    MD_sine_wave(out, 64, 1.0, 0.0, 44100.0);
    for (unsigned i = 0; i < 64; i++) {
        if (fabs(out[i]) > 1e-15) return 0;
    }
    return 1;
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

/** Zero amplitude produces all-zero output. */
static int test_white_noise_zero_amplitude(void)
{
    double out[128];
    MD_white_noise(out, 128, 0.0, 42);
    for (unsigned i = 0; i < 128; i++) {
        if (fabs(out[i]) > 1e-15) return 0;
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

/** Zero amplitude produces all-zero output. */
static int test_square_zero_amplitude(void)
{
    double out[64];
    MD_square_wave(out, 64, 0.0, 1000.0, 16000.0);
    for (unsigned i = 0; i < 64; i++) {
        if (fabs(out[i]) > 1e-15) return 0;
    }
    return 1;
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

/** Zero amplitude produces all-zero output. */
static int test_sawtooth_zero_amplitude(void)
{
    double out[64];
    MD_sawtooth_wave(out, 64, 0.0, 1000.0, 16000.0);
    for (unsigned i = 0; i < 64; i++) {
        if (fabs(out[i]) > 1e-15) return 0;
    }
    return 1;
}

/* -----------------------------------------------------------------------
 * MD_shepard_tone
 * -----------------------------------------------------------------------*/

/** Output has non-zero energy. */
static int test_shepard_nonzero_energy(void)
{
    unsigned N = 44100;
    double *buf = malloc(N * sizeof(double));
    MD_shepard_tone(buf, N, 0.8, 440.0, 44100.0, 0.5, 8);
    double energy = 0.0;
    for (unsigned i = 0; i < N; i++) energy += buf[i] * buf[i];
    free(buf);
    return energy > 0.0;
}

/** Peak absolute amplitude matches the requested amplitude. */
static int test_shepard_peak_amplitude(void)
{
    unsigned N = 44100;
    double *buf = malloc(N * sizeof(double));
    MD_shepard_tone(buf, N, 0.8, 440.0, 44100.0, 0.5, 8);
    double peak = 0.0;
    for (unsigned i = 0; i < N; i++) {
        double a = fabs(buf[i]);
        if (a > peak) peak = a;
    }
    free(buf);
    return approx_equal(peak, 0.8, 0.01);
}

/** Static Shepard chord (rate=0) should have octave-spaced spectral peaks.
 *  Use 4 octaves centred on 1000 Hz: expect peaks near 250, 500, 1000, 2000 Hz. */
static int test_shepard_static_octave_peaks(void)
{
    /* Use N = fs so that bin k = k Hz exactly (no spectral leakage). */
    unsigned N = 8192;
    double sr = 8192.0;
    double base = 1024.0;
    unsigned num_oct = 4;

    double *buf = malloc(N * sizeof(double));
    double *mag = malloc((N / 2 + 1) * sizeof(double));

    MD_shepard_tone(buf, N, 1.0, base, sr, 0.0, num_oct);
    MD_magnitude_spectrum(buf, N, mag);

    /* Expected peaks: base * 2^(k - center) for k = 0..3, center = 1.5
     * k=0: 1024 * 2^(-1.5) = 362.0
     * k=1: 1024 * 2^(-0.5) = 724.1
     * k=2: 1024 * 2^( 0.5) = 1448.2
     * k=3: 1024 * 2^( 1.5) = 2896.3 */
    double expected[] = {362.0, 724.1, 1448.2, 2896.3};
    int ok = 1;
    for (int e = 0; e < 4; e++) {
        unsigned bin = (unsigned)(expected[e] + 0.5); /* bin k ≈ freq Hz */
        /* Check that this bin is a local maximum above threshold */
        double thresh = mag[1]; /* arbitrary low bin as baseline */
        for (unsigned b = 1; b < N / 2 + 1; b++)
            if (mag[b] < thresh) thresh = mag[b];
        /* The peak bin should be notably above the noise floor */
        int found = 0;
        for (int delta = -3; delta <= 3; delta++) {
            int b = (int)bin + delta;
            if (b < 0 || b >= (int)(N / 2 + 1)) continue;
            if (mag[b] > thresh * 10.0) found = 1;
        }
        if (!found) ok = 0;
    }
    free(mag);
    free(buf);
    return ok;
}

/** Rising and falling Shepard tones should differ. */
static int test_shepard_rising_vs_falling(void)
{
    unsigned N = 22050; /* 0.5 s at 44100 Hz */
    double *rising  = malloc(N * sizeof(double));
    double *falling = malloc(N * sizeof(double));

    MD_shepard_tone(rising,  N, 0.8, 440.0, 44100.0,  0.5, 8);
    MD_shepard_tone(falling, N, 0.8, 440.0, 44100.0, -0.5, 8);

    /* They should not be identical */
    int differ = 0;
    for (unsigned i = 0; i < N; i++) {
        if (fabs(rising[i] - falling[i]) > 1e-10) { differ = 1; break; }
    }
    free(rising);
    free(falling);
    return differ;
}

/** Gaussian spectral envelope: the layer nearest the centre should be
 *  louder in the spectrum than layers at the edge.
 *
 *  Use odd num_octaves so that one layer sits exactly at base_freq.
 *  With N = fs = 8192, bin k = k Hz exactly. */
static int test_shepard_gaussian_envelope(void)
{
    unsigned N = 8192;
    double sr = 8192.0;
    double base = 512.0;
    unsigned num_oct = 7;  /* centre = 3.0 → layer 3 at exactly 512 Hz */

    double *buf = malloc(N * sizeof(double));
    double *mag = malloc((N / 2 + 1) * sizeof(double));

    MD_shepard_tone(buf, N, 1.0, base, sr, 0.0, num_oct);
    MD_magnitude_spectrum(buf, N, mag);

    /* Layer 3: d=0, freq = 512 Hz (bin 512), Gaussian = 1.0
     * Layer 0: d=-3, freq = 64 Hz (bin 64), Gaussian ≈ 0.23
     * The peak at the centre should be substantially larger than the edge. */
    double centre_peak = 0.0;
    for (int d = -2; d <= 2; d++) {
        int b = 512 + d;
        if (b >= 0 && b < (int)(N / 2 + 1) && mag[b] > centre_peak)
            centre_peak = mag[b];
    }
    double edge_peak = 0.0;
    for (int d = -2; d <= 2; d++) {
        int b = 64 + d;
        if (b >= 0 && b < (int)(N / 2 + 1) && mag[b] > edge_peak)
            edge_peak = mag[b];
    }

    int ok = (centre_peak > edge_peak * 2.0);
    free(mag);
    free(buf);
    return ok;
}

/** Zero amplitude produces all-zero output. */
static int test_shepard_zero_amplitude(void)
{
    unsigned N = 8192;
    double *buf = malloc(N * sizeof(double));
    MD_shepard_tone(buf, N, 0.0, 440.0, 44100.0, 0.5, 8);
    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        if (fabs(buf[i]) > 1e-15) { ok = 0; break; }
    }
    free(buf);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_generators_tests(void)
{
    printf("\n--- MD_sine_wave ---\n");
    RUN_TEST(test_sine_starts_at_zero);
    RUN_TEST(test_sine_quarter_period);
    RUN_TEST(test_sine_full_period);
    RUN_TEST(test_sine_spectrum_peak);
    RUN_TEST(test_sine_amplitude_negative);
    RUN_TEST(test_sine_zero_amplitude);
    RUN_TEST(test_sine_zero_frequency);

    printf("\n--- MD_white_noise ---\n");
    RUN_TEST(test_white_noise_mean_near_zero);
    RUN_TEST(test_white_noise_stddev_matches_amplitude);
    RUN_TEST(test_white_noise_reproducible);
    RUN_TEST(test_white_noise_different_seeds);
    RUN_TEST(test_white_noise_flat_spectrum);
    RUN_TEST(test_white_noise_odd_length);
    RUN_TEST(test_white_noise_zero_amplitude);

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
    RUN_TEST(test_square_zero_amplitude);

    printf("\n--- MD_sawtooth_wave ---\n");
    RUN_TEST(test_sawtooth_midpoint);
    RUN_TEST(test_sawtooth_spectrum_harmonics);
    RUN_TEST(test_sawtooth_amplitude_negative);
    RUN_TEST(test_sawtooth_zero_amplitude);

    printf("\n--- MD_shepard_tone ---\n");
    RUN_TEST(test_shepard_nonzero_energy);
    RUN_TEST(test_shepard_peak_amplitude);
    RUN_TEST(test_shepard_static_octave_peaks);
    RUN_TEST(test_shepard_rising_vs_falling);
    RUN_TEST(test_shepard_gaussian_envelope);
    RUN_TEST(test_shepard_zero_amplitude);
}
