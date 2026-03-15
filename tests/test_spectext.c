/**
 * @file test_spectext.c
 * @brief Tests for MD_spectrogram_text.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * MD_spectrogram_text
 * -----------------------------------------------------------------------*/

/** Output length matches expected col_samples * grid_cols. */
static int test_spectext_output_length(void)
{
    const char *text = "AB";
    double sr = 16000.0;
    double dur = 1.0;
    unsigned len = (unsigned)strlen(text);
    unsigned grid_cols = len * 8 - 3;                     /* 13 */
    unsigned col_samples = (unsigned)(dur / (double)grid_cols * sr);
    unsigned expected = col_samples * grid_cols;

    double *buf = malloc(expected * sizeof(double));
    unsigned n = MD_spectrogram_text(buf, expected, text,
                                     200.0, 7500.0, dur, sr);
    int ok = (n == expected);
    free(buf);
    return ok;
}

/** Single character "A": verify non-zero energy in the output. */
static int test_spectext_nonzero_energy(void)
{
    unsigned buflen = 32000;
    double *buf = malloc(buflen * sizeof(double));
    unsigned n = MD_spectrogram_text(buf, buflen, "A",
                                     200.0, 7500.0, 2.0, 16000.0);
    double energy = 0.0;
    for (unsigned i = 0; i < n; i++) energy += buf[i] * buf[i];
    free(buf);
    return energy > 0.0;
}

/** Frequency check: "I" has all 7 rows on in its center columns; verify
 *  the spectral peak is within the expected frequency range. */
static int test_spectext_frequency_peak(void)
{
    double sr = 16000.0;
    double freq_lo = 1000.0;
    double freq_hi = 7000.0;
    unsigned buflen = 64000;
    double *buf = malloc(buflen * sizeof(double));
    /* Use "I" - has a vertical bar at columns 1-3 (rows 0 through 6 all on) */
    unsigned n = MD_spectrogram_text(buf, buflen, "I",
                                     freq_lo, freq_hi, 2.0, sr);

    /* Compute magnitude spectrum of a chunk from the middle */
    unsigned fft_n = 4096;
    if (n < fft_n) { free(buf); return 0; }
    unsigned offset = (n - fft_n) / 2;

    double *mag = malloc((fft_n / 2 + 1) * sizeof(double));
    MD_magnitude_spectrum(buf + offset, fft_n, mag);

    /* Row 0 = freq_hi, row 6 = freq_lo.  "I" has all 7 rows on in its
     * center column.  Find the overall peak — should be at one of the
     * 7 row frequencies. */
    unsigned peak_bin = 1;
    for (unsigned k = 2; k < fft_n / 2 + 1; k++) {
        if (mag[k] > mag[peak_bin]) peak_bin = k;
    }
    double peak_freq = (double)peak_bin * sr / (double)fft_n;

    /* Check that the peak is between freq_lo and freq_hi */
    int ok = (peak_freq >= freq_lo - 100.0 && peak_freq <= freq_hi + 100.0);

    free(mag);
    free(buf);
    return ok;
}

/** Normalization: verify peak absolute value is approximately 0.9. */
static int test_spectext_normalization(void)
{
    unsigned buflen = 64000;
    double *buf = malloc(buflen * sizeof(double));
    unsigned n = MD_spectrogram_text(buf, buflen, "HELLO",
                                     200.0, 7500.0, 2.0, 16000.0);

    double peak = 0.0;
    for (unsigned i = 0; i < n; i++) {
        double a = fabs(buf[i]);
        if (a > peak) peak = a;
    }
    free(buf);
    return approx_equal(peak, 0.9, 0.01);
}

/** Space character produces silence (all zeros). */
static int test_spectext_space_silence(void)
{
    unsigned buflen = 32000;
    double *buf = malloc(buflen * sizeof(double));
    unsigned n = MD_spectrogram_text(buf, buflen, " ",
                                     200.0, 7500.0, 2.0, 16000.0);

    int all_zero = 1;
    for (unsigned i = 0; i < n; i++) {
        if (buf[i] != 0.0) { all_zero = 0; break; }
    }
    free(buf);
    return all_zero;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_spectext_tests(void)
{
    printf("\n--- MD_spectrogram_text ---\n");
    RUN_TEST(test_spectext_output_length);
    RUN_TEST(test_spectext_nonzero_energy);
    RUN_TEST(test_spectext_frequency_peak);
    RUN_TEST(test_spectext_normalization);
    RUN_TEST(test_spectext_space_silence);
}
