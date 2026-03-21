/**
 * @file test_vad.c
 * @brief Tests for the VAD module.
 */

#include <math.h>
#include "minidsp.h"
#include "fileio.h"
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
 * Feature directionality tests
 * -----------------------------------------------------------------------*/

/**
 * Verify that all features are higher for a speech-like signal (tonal)
 * than for a noise-like signal (flat spectrum).  This catches
 * directionality bugs where high entropy/flatness (noise traits)
 * would incorrectly push the score toward "speech".
 */
static int test_vad_feature_directionality(void)
{
    MD_vad_params p;
    MD_vad_default_params(&p);
    p.adaptation_rate = 0.5; /* fast adaptation so min/max converges */

    MD_vad_state st;
    MD_vad_init(&st, &p);

    unsigned N = 256;
    double tone[256];
    double noise[256];

    /* Speech-like: a 1 kHz sine (structured spectrum) */
    MD_sine_wave(tone, N, 0.5, 1000.0, 16000.0);

    /* Noise-like: white noise (flat spectrum) */
    MD_white_noise(noise, N, 0.5, 42);

    /* Calibrate on alternating signals so min/max brackets both */
    for (int i = 0; i < 100; i++) {
        MD_vad_calibrate(&st, tone, N, 16000.0);
        MD_vad_calibrate(&st, noise, N, 16000.0);
    }

    /* Process one tone frame, capture features */
    double feat_tone[MD_VAD_NUM_FEATURES];
    double score_tone = 0.0;
    MD_vad_process_frame(&st, tone, N, 16000.0, &score_tone, feat_tone);

    /* Process one noise frame, capture features */
    double feat_noise[MD_VAD_NUM_FEATURES];
    double score_noise = 0.0;
    MD_vad_process_frame(&st, noise, N, 16000.0, &score_noise, feat_noise);

    /* Spectral entropy: should be HIGHER for tone (speech) than noise.
     * Before the inversion fix, this would fail. */
    if (feat_tone[MD_VAD_FEAT_SPECTRAL_ENTROPY] <=
        feat_noise[MD_VAD_FEAT_SPECTRAL_ENTROPY])
        return 0;

    /* Spectral flatness: should be HIGHER for tone (speech) than noise. */
    if (feat_tone[MD_VAD_FEAT_SPECTRAL_FLATNESS] <=
        feat_noise[MD_VAD_FEAT_SPECTRAL_FLATNESS])
        return 0;

    /* Overall score should be higher for the speech-like signal */
    if (score_tone <= score_noise) return 0;

    return 1;
}

/* -----------------------------------------------------------------------
 * Real speech test (sa2.wav)
 * -----------------------------------------------------------------------*/

/**
 * Run VAD on a real speech recording and verify it matches known
 * speech boundaries.  sa2.wav has speech from ~0.108s to ~2.167s
 * with silence before and after.
 */
static int test_vad_real_speech(void)
{
    float *audio_f = NULL;
    size_t num_samples = 0;
    unsigned sample_rate = 0;

    if (FIO_read_audio("../samples/sa2.wav", &audio_f, &num_samples,
                       &sample_rate, 1) != 0) {
        printf("(skip: sa2.wav not found) ");
        return 1;
    }

    /* Convert float to double for VAD API */
    double *audio = malloc(num_samples * sizeof(double));
    if (!audio) { free(audio_f); return 0; }
    for (size_t i = 0; i < num_samples; i++)
        audio[i] = (double)audio_f[i];
    free(audio_f);

    unsigned N = 256;   /* 16ms frame at 16 kHz */
    unsigned hop = 256;
    double sr = (double)sample_rate;

    MD_vad_params p;
    MD_vad_default_params(&p);
    p.onset_frames = 3;
    p.hangover_frames = 3; /* short hangover so trailing silence is detected */
    MD_vad_state st;
    MD_vad_init(&st, &p);

    /* Calibrate on silence before speech onset (first ~0.1s) */
    double speech_start_s = 0.108;
    unsigned calib_end = (unsigned)(speech_start_s * sr);
    for (unsigned off = 0; off + N <= calib_end; off += hop)
        MD_vad_calibrate(&st, audio + off, N, sr);

    /* Process all frames, record decisions */
    unsigned num_frames = (num_samples >= N) ? (num_samples - N) / hop + 1 : 0;
    int *decisions = malloc(num_frames * sizeof(int));
    if (!decisions) { free(audio); return 0; }

    for (unsigned f = 0; f < num_frames; f++) {
        unsigned off = f * hop;
        decisions[f] = MD_vad_process_frame(&st, audio + off, N, sr,
                                            NULL, NULL);
    }

    int ok = 1;

    /* 1. Pre-speech silence: frame at 0.05s should be inactive */
    unsigned pre_frame = (unsigned)(0.05 * sr / hop);
    if (pre_frame < num_frames && decisions[pre_frame] != 0) ok = 0;

    /* 2. Mid-speech: frame at 1.0s should be active */
    unsigned mid_frame = (unsigned)(1.0 * sr / hop);
    if (mid_frame < num_frames && decisions[mid_frame] != 1) ok = 0;

    /* 3. Trailing silence: last frame should be inactive */
    if (num_frames > 0 && decisions[num_frames - 1] != 0) ok = 0;

    free(decisions);
    free(audio);
    return ok;
}

/**
 * Run VAD on a second real speech recording.  berp_01_1_0001.wav has
 * speech from ~0.446s to ~3.235s with silence before and after.
 */
static int test_vad_real_speech_berp(void)
{
    float *audio_f = NULL;
    size_t num_samples = 0;
    unsigned sample_rate = 0;

    if (FIO_read_audio("../samples/berp_01_1_0001.wav", &audio_f,
                       &num_samples, &sample_rate, 1) != 0) {
        printf("(skip: berp not found) ");
        return 1;
    }

    double *audio = malloc(num_samples * sizeof(double));
    if (!audio) { free(audio_f); return 0; }
    for (size_t i = 0; i < num_samples; i++)
        audio[i] = (double)audio_f[i];
    free(audio_f);

    unsigned N = 256;
    unsigned hop = 256;
    double sr = (double)sample_rate;

    MD_vad_params p;
    MD_vad_default_params(&p);
    p.onset_frames = 3;
    p.hangover_frames = 3;
    p.threshold = 0.6; /* above noise baseline (~0.5 after calibration) */
    MD_vad_state st;
    MD_vad_init(&st, &p);

    /* Calibrate on silence before speech onset (first ~0.4s) */
    unsigned calib_end = (unsigned)(0.4 * sr);
    for (unsigned off = 0; off + N <= calib_end; off += hop)
        MD_vad_calibrate(&st, audio + off, N, sr);

    unsigned num_frames = (num_samples >= N) ? (num_samples - N) / hop + 1 : 0;
    int *decisions = malloc(num_frames * sizeof(int));
    if (!decisions) { free(audio); return 0; }

    for (unsigned f = 0; f < num_frames; f++) {
        unsigned off = f * hop;
        decisions[f] = MD_vad_process_frame(&st, audio + off, N, sr,
                                            NULL, NULL);
    }

    int ok = 1;

    /* 1. Pre-speech silence: frame at 0.35s should be inactive.
     *    (Avoid the first ~0.25s where adaptive normalization settles.) */
    unsigned pre_frame = (unsigned)(0.35 * sr / hop);
    if (pre_frame < num_frames && decisions[pre_frame] != 0) ok = 0;

    /* 2. Mid-speech: frame at 1.8s should be active */
    unsigned mid_frame = (unsigned)(1.8 * sr / hop);
    if (mid_frame < num_frames && decisions[mid_frame] != 1) ok = 0;

    /* 3. Trailing silence: last frame should be inactive */
    if (num_frames > 0 && decisions[num_frames - 1] != 0) ok = 0;

    free(decisions);
    free(audio);
    return ok;
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
    RUN_TEST(test_vad_feature_directionality);
    RUN_TEST(test_vad_real_speech);
    RUN_TEST(test_vad_real_speech_berp);
}
