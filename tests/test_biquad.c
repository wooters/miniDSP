/**
 * @file test_biquad.c
 * @brief Tests for BiQuad filters.
 */

#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "biquad.h"
#include "test_helpers.h"

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

/** BiQuad_new should return NULL for an invalid filter type. */
static int test_biquad_invalid_type(void)
{
    biquad *b = BiQuad_new(999, 0.0, 1000.0, 44100.0, 1.0);
    return (b == NULL);
}

/** BiQuad_new should return NULL for freq=0 (sin(omega)==0). */
static int test_biquad_freq_zero(void)
{
    biquad *b = BiQuad_new(LPF, 0.0, 0.0, 44100.0, 1.0);
    return (b == NULL);
}

/** BiQuad_new should return NULL for freq=Nyquist (sin(omega)==0). */
static int test_biquad_freq_nyquist(void)
{
    biquad *b = BiQuad_new(LPF, 0.0, 22050.0, 44100.0, 1.0);
    return (b == NULL);
}

/** Frequencies just inside valid range should succeed. */
static int test_biquad_freq_near_boundaries(void)
{
    biquad *b1 = BiQuad_new(LPF, 0.0, 1.0, 44100.0, 1.0);
    biquad *b2 = BiQuad_new(LPF, 0.0, 22049.0, 44100.0, 1.0);
    int ok = (b1 != NULL) && (b2 != NULL);
    free(b2);
    free(b1);
    return ok;
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

void run_biquad_tests(void)
{
    printf("\n--- Biquad filters ---\n");
    RUN_TEST(test_biquad_invalid_type);
    RUN_TEST(test_biquad_freq_zero);
    RUN_TEST(test_biquad_freq_nyquist);
    RUN_TEST(test_biquad_freq_near_boundaries);
    RUN_TEST(test_biquad_lpf);
    RUN_TEST(test_biquad_hpf);
    RUN_TEST(test_biquad_bpf);
    RUN_TEST(test_biquad_notch);
    RUN_TEST(test_biquad_peq_boost);
    RUN_TEST(test_biquad_low_shelf);
    RUN_TEST(test_biquad_high_shelf);
    RUN_TEST(test_biquad_hpf_dc_rejection);
}
