/**
 * @file test_effects.c
 * @brief Tests for simple audio effects (delay echo, tremolo, comb reverb).
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

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

/** Tremolo at full depth (1.0) should reach zero gain. */
static int test_tremolo_full_depth(void)
{
    const unsigned N = 8000;
    double in[N], out[N];
    for (unsigned i = 0; i < N; i++) in[i] = 1.0;

    MD_tremolo(in, out, N, 5.0, 1.0, 8000.0);

    /* With depth=1.0, gain oscillates between 0.0 and 1.0.
     * Find the minimum output value — it should reach ~0. */
    double min_val = out[0];
    for (unsigned i = 1; i < N; i++) {
        if (out[i] < min_val) min_val = out[i];
    }
    return approx_equal(min_val, 0.0, 0.01);
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

/** Comb reverb with feedback=0: impulse produces input + one delayed copy, no further repeats. */
static int test_comb_reverb_zero_feedback(void)
{
    const unsigned N = 20;
    double in[N], out[N];
    memset(in, 0, sizeof(in));
    in[0] = 1.0;

    MD_comb_reverb(in, out, N, 4, 0.0, 0.0, 1.0);

    /* With feedback=0: y_comb = x + 0*delayed = x.
     * So wet output is just the input. Only impulse at index 0. */
    int ok = 1;
    ok &= approx_equal(out[0], 1.0, 1e-12);
    for (unsigned i = 1; i < N; i++) {
        ok &= approx_equal(out[i], 0.0, 1e-12);
    }
    return ok;
}

void run_effects_tests(void)
{
    printf("\n--- Simple effects ---\n");
    RUN_TEST(test_delay_echo_dry_passthrough);
    RUN_TEST(test_delay_echo_impulse_decay);
    RUN_TEST(test_delay_echo_inplace);
    RUN_TEST(test_tremolo_depth_zero_passthrough);
    RUN_TEST(test_tremolo_gain_bounds);
    RUN_TEST(test_tremolo_full_depth);
    RUN_TEST(test_comb_reverb_dry_passthrough);
    RUN_TEST(test_comb_reverb_impulse_decay);
    RUN_TEST(test_comb_reverb_inplace);
    RUN_TEST(test_comb_reverb_zero_feedback);
}
