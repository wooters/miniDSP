/**
 * @file test_windows.c
 * @brief Tests for window generation functions.
 */

#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

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

/** Hann window of length 2: both endpoints should be zero. */
static int test_hann_length_two(void)
{
    double win[2];
    MD_Gen_Hann_Win(win, 2);
    return approx_equal(win[0], 0.0, 1e-15)
        && approx_equal(win[1], 0.0, 1e-15);
}

/** Hamming window of length 2: both endpoints equal the minimum (~0.08). */
static int test_hamming_length_two(void)
{
    double win[2];
    MD_Gen_Hamming_Win(win, 2);
    /* Hamming: w[n] = 0.54 - 0.46*cos(2*pi*n/(N-1)). For N=2: n=0 → 0.08, n=1 → 0.08 */
    int ok = approx_equal(win[0], 0.08, 1e-10);
    ok &= approx_equal(win[1], 0.08, 1e-10);
    return ok;
}

/** Blackman window of length 2: both endpoints near zero. */
static int test_blackman_length_two(void)
{
    double win[2];
    MD_Gen_Blackman_Win(win, 2);
    /* Blackman: w[n] = 0.42 - 0.5*cos(2pi*n/(N-1)) + 0.08*cos(4pi*n/(N-1)).
     * For N=2: n=0 → 0.42 - 0.5 + 0.08 = 0.0, n=1 → same */
    return approx_equal(win[0], 0.0, 1e-10)
        && approx_equal(win[1], 0.0, 1e-10);
}

/* -----------------------------------------------------------------------
 * Runner
 * -----------------------------------------------------------------------*/

void run_windows_tests(void)
{
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
    RUN_TEST(test_hann_length_two);
    RUN_TEST(test_hamming_length_two);
    RUN_TEST(test_blackman_length_two);
}
