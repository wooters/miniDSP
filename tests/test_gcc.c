/**
 * @file test_gcc.c
 * @brief Tests for GCC-PHAT delay estimation.
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

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

    int delay = MD_get_delay(siga, sigb, N, NULL, 50, PHAT);
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

    int delay = MD_get_delay(siga, sigb, N, NULL, 50, PHAT);
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
    int delay = MD_get_delay(siga, siga, N, NULL, 50, PHAT);
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

    int delay = MD_get_delay(siga, sigb, N, NULL, 50, SIMP);
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
        manual_results[k] = MD_get_delay(w_ref, w_sig, N, NULL, 50, PHAT);
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
        int d = MD_get_delay(siga, sigb, N, NULL, 50, PHAT);
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
        int d = MD_get_delay(siga, sigb, N, NULL, 50, PHAT);
        ok &= (d == true_delay);
        free(sigb);
        free(siga);
    }

    return ok;
}

/** Identical signals should produce zero delay estimate. */
static int test_gcc_identical_signals(void)
{
    unsigned N = 4096;
    double *sig = calloc(N, sizeof(double));
    for (unsigned i = 0; i < N; i++)
        sig[i] = sin(2.0 * M_PI * (double)i / 64.0);
    int d = MD_get_delay(sig, sig, N, NULL, 50, PHAT);
    free(sig);
    return (d == 0);
}

/** Two silent signals should produce zero cross-correlation. */
static int test_gcc_silence(void)
{
    unsigned N = 1024;
    double *siga = calloc(N, sizeof(double));
    double *sigb = calloc(N, sizeof(double));
    double *lagvals = malloc(N * sizeof(double));

    MD_gcc(siga, sigb, N, lagvals, PHAT);

    int ok = 1;
    for (unsigned i = 0; i < N; i++) {
        ok &= approx_equal(lagvals[i], 0.0, 1e-15);
    }

    free(lagvals);
    free(sigb);
    free(siga);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public runner
 * -----------------------------------------------------------------------*/

void run_gcc_tests(void)
{
    printf("\n--- GCC-PHAT delay estimation ---\n");
    RUN_TEST(test_gcc_phat_positive_delay);
    RUN_TEST(test_gcc_phat_negative_delay);
    RUN_TEST(test_gcc_phat_zero_delay);
    RUN_TEST(test_gcc_simp_delay);
    RUN_TEST(test_gcc_returns_entropy);
    RUN_TEST(test_gcc_multiple_delays);
    RUN_TEST(test_gcc_lagvals_structure);
    RUN_TEST(test_gcc_different_lengths);
    RUN_TEST(test_gcc_identical_signals);
    RUN_TEST(test_gcc_silence);
}
