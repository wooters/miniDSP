/**
 * @file test_vad.c
 * @brief Tests for the VAD module.
 */

#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * MD_vad_default_params tests
 * -----------------------------------------------------------------------*/

static int test_vad_default_params(void)
{
    MD_vad_params p;
    MD_vad_default_params(&p);

    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++) {
        if (!approx_equal(p.weights[i], 0.2, 1e-12)) return 0;
    }
    if (!approx_equal(p.threshold, 0.5, 1e-12)) return 0;
    if (p.onset_frames != 3) return 0;
    if (p.hangover_frames != 15) return 0;
    if (!approx_equal(p.adaptation_rate, 0.01, 1e-12)) return 0;
    if (!approx_equal(p.band_low_hz, 300.0, 1e-12)) return 0;
    if (!approx_equal(p.band_high_hz, 3400.0, 1e-12)) return 0;

    return 1;
}

/* -----------------------------------------------------------------------
 * MD_vad_calibrate tests
 * -----------------------------------------------------------------------*/

static int test_vad_calibrate(void)
{
    MD_vad_state st;
    MD_vad_init(&st, NULL);

    unsigned N = 256;
    double silence[256] = {0};

    /* Calibrate on 20 frames of silence */
    for (int i = 0; i < 20; i++)
        MD_vad_calibrate(&st, silence, N, 16000.0);

    /* After calibration, min/max should be close to silence values */
    if (st.frames_processed != 20) return 0;

    /* Silence energy is 0, min/max should both be near 0 for energy */
    if (!approx_equal(st.feat_min[MD_VAD_FEAT_ENERGY],
                      st.feat_max[MD_VAD_FEAT_ENERGY], 1e-6)) return 0;

    /* Decision should still be 0 (calibrate doesn't run state machine) */
    if (st.current_decision != 0) return 0;
    if (st.onset_counter != 0) return 0;

    return 1;
}

/* -----------------------------------------------------------------------
 * MD_vad_process_frame tests
 * -----------------------------------------------------------------------*/

static int test_vad_silence(void)
{
    MD_vad_state st;
    MD_vad_init(&st, NULL);

    unsigned N = 256;
    double silence[256] = {0};

    /* Process 30 frames of silence */
    double score = 0.0;
    int decision = 0;
    for (int i = 0; i < 30; i++) {
        decision = MD_vad_process_frame(&st, silence, N, 16000.0,
                                        &score, NULL);
    }

    /* Should remain silence */
    if (decision != 0) return 0;

    return 1;
}

static int test_vad_speech(void)
{
    MD_vad_state st;
    MD_vad_params p;
    MD_vad_default_params(&p);
    p.onset_frames = 3;
    p.threshold = 0.3;
    MD_vad_init(&st, &p);

    unsigned N = 256;
    double silence[256] = {0};
    double speech[256];

    /* Generate speech-like signal: sine at 1000 Hz */
    MD_sine_wave(speech, N, 0.5, 1000.0, 16000.0);

    /* Calibrate on silence first */
    for (int i = 0; i < 20; i++)
        MD_vad_calibrate(&st, silence, N, 16000.0);

    /* Process enough speech frames to trigger onset */
    int decision = 0;
    for (int i = 0; i < 30; i++) {
        decision = MD_vad_process_frame(&st, speech, N, 16000.0,
                                        NULL, NULL);
    }

    /* Should detect speech */
    if (decision != 1) return 0;

    return 1;
}

static int test_vad_hangover(void)
{
    MD_vad_state st;
    MD_vad_params p;
    MD_vad_default_params(&p);
    p.onset_frames = 2;
    p.hangover_frames = 5;
    p.threshold = 0.3;
    MD_vad_init(&st, &p);

    unsigned N = 256;
    double silence[256] = {0};
    double speech[256];
    MD_sine_wave(speech, N, 0.5, 1000.0, 16000.0);

    /* Calibrate on silence */
    for (int i = 0; i < 20; i++)
        MD_vad_calibrate(&st, silence, N, 16000.0);

    /* Speech frames to trigger detection */
    int decision = 0;
    for (int i = 0; i < 20; i++)
        decision = MD_vad_process_frame(&st, speech, N, 16000.0,
                                        NULL, NULL);
    if (decision != 1) return 0;

    /* Now send silence — decision should hold for hangover_frames */
    for (unsigned i = 0; i < p.hangover_frames; i++) {
        decision = MD_vad_process_frame(&st, silence, N, 16000.0,
                                        NULL, NULL);
        /* During hangover, decision should still be 1
         * (it drops to 0 when counter hits 0, which is at the end) */
    }

    /* After hangover_frames of silence, decision should drop to 0 */
    /* Give a few more frames to ensure it transitions */
    for (int i = 0; i < 5; i++)
        decision = MD_vad_process_frame(&st, silence, N, 16000.0,
                                        NULL, NULL);

    if (decision != 0) return 0;

    return 1;
}

static int test_vad_onset(void)
{
    MD_vad_state st;
    MD_vad_params p;
    MD_vad_default_params(&p);
    p.onset_frames = 10;
    p.threshold = 0.3;
    MD_vad_init(&st, &p);

    unsigned N = 256;
    double silence[256] = {0};
    double speech[256];
    MD_sine_wave(speech, N, 0.5, 1000.0, 16000.0);

    /* Calibrate on silence */
    for (int i = 0; i < 20; i++)
        MD_vad_calibrate(&st, silence, N, 16000.0);

    /* Send fewer than onset_frames speech frames */
    int decision = 0;
    for (int i = 0; i < 5; i++)
        decision = MD_vad_process_frame(&st, speech, N, 16000.0,
                                        NULL, NULL);

    /* Should NOT have triggered speech (need 10 consecutive frames) */
    if (decision != 0) return 0;

    return 1;
}

static int test_vad_features_out(void)
{
    MD_vad_state st;
    MD_vad_init(&st, NULL);

    unsigned N = 256;
    double speech[256];
    MD_sine_wave(speech, N, 0.5, 1000.0, 16000.0);

    /* Process several frames so normalization stabilizes */
    double features[MD_VAD_NUM_FEATURES];
    for (int i = 0; i < 20; i++)
        MD_vad_process_frame(&st, speech, N, 16000.0, NULL, features);

    /* All features should be in [0.0, 1.0] */
    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++) {
        if (features[i] < 0.0 || features[i] > 1.0)
            return 0;
    }

    return 1;
}

static int test_vad_custom_weights(void)
{
    MD_vad_state st;
    MD_vad_params p;
    MD_vad_default_params(&p);

    /* Set only energy weight to 1.0, rest to 0 */
    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++)
        p.weights[i] = 0.0;
    p.weights[MD_VAD_FEAT_ENERGY] = 1.0;
    MD_vad_init(&st, &p);

    unsigned N = 256;
    double silence[256] = {0};
    double speech[256];
    MD_sine_wave(speech, N, 0.5, 1000.0, 16000.0);

    /* Calibrate on silence */
    for (int i = 0; i < 10; i++)
        MD_vad_calibrate(&st, silence, N, 16000.0);

    /* Process speech frame and check score tracks energy */
    double score_energy = 0.0;
    double features[MD_VAD_NUM_FEATURES];
    MD_vad_process_frame(&st, speech, N, 16000.0, &score_energy, features);

    /* Score should equal the energy feature since only energy has weight */
    if (!approx_equal(score_energy, features[MD_VAD_FEAT_ENERGY], 1e-12))
        return 0;

    return 1;
}

/* -----------------------------------------------------------------------
 * Test runner
 * -----------------------------------------------------------------------*/

void run_vad_tests(void)
{
    printf("\n--- VAD ---\n");
    RUN_TEST(test_vad_default_params);
    RUN_TEST(test_vad_calibrate);
    RUN_TEST(test_vad_silence);
    RUN_TEST(test_vad_speech);
    RUN_TEST(test_vad_hangover);
    RUN_TEST(test_vad_onset);
    RUN_TEST(test_vad_features_out);
    RUN_TEST(test_vad_custom_weights);
}
