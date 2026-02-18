/**
 * @file gen_audio_samples.c
 * @brief Docs-only utility: generate WAV audio samples for docs guides.
 *
 * This program generates WAV files into guides/audio/ so that the Doxygen
 * documentation can embed playable audio previews for signal generators and
 * simple effects.
 *
 * It is NOT a user-facing example — it is invoked automatically by `make docs`.
 *
 * Build (from repo root):
 *   make -C examples gen_audio_samples
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "fileio.h"

#define SAMPLE_RATE 44100
#define DURATION    2.0

static double clamp_unit(double x)
{
    if (x > 1.0) return 1.0;
    if (x < -1.0) return -1.0;
    return x;
}

static int write_wav(const char *path, const double *buf, unsigned N)
{
    float *fbuf = malloc(N * sizeof(float));
    if (!fbuf) {
        fprintf(stderr, "allocation failed for %s\n", path);
        return -1;
    }
    for (unsigned i = 0; i < N; i++)
        fbuf[i] = (float)clamp_unit(buf[i]);

    int rc = FIO_write_wav(path, fbuf, N, SAMPLE_RATE);
    free(fbuf);
    if (rc != 0)
        fprintf(stderr, "failed to write %s\n", path);
    else
        printf("  %s\n", path);
    return rc;
}

int main(void)
{
    const unsigned N = (unsigned)(SAMPLE_RATE * DURATION);
    double *buf = malloc(N * sizeof(double));
    double *fx  = malloc(N * sizeof(double));
    if (!buf || !fx) {
        fprintf(stderr, "allocation failed\n");
        free(fx);
        free(buf);
        return 1;
    }

    printf("Generating audio samples (%u samples, %.0f Hz, %.1fs):\n",
           N, (double)SAMPLE_RATE, DURATION);

    /* Sine wave — 440 Hz */
    MD_sine_wave(buf, N, 0.8, 440.0, SAMPLE_RATE);
    write_wav("guides/audio/sine_440hz.wav", buf, N);

    /* Square wave — 440 Hz */
    MD_square_wave(buf, N, 0.8, 440.0, SAMPLE_RATE);
    write_wav("guides/audio/square_440hz.wav", buf, N);

    /* Sawtooth wave — 440 Hz */
    MD_sawtooth_wave(buf, N, 0.8, 440.0, SAMPLE_RATE);
    write_wav("guides/audio/sawtooth_440hz.wav", buf, N);

    /* Linear chirp — 20 Hz to 4000 Hz */
    MD_chirp_linear(buf, N, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_wav("guides/audio/chirp_linear.wav", buf, N);

    /* Logarithmic chirp — 20 Hz to 4000 Hz */
    MD_chirp_log(buf, N, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_wav("guides/audio/chirp_log.wav", buf, N);

    /* White noise — sigma 0.25, seed 42 */
    MD_white_noise(buf, N, 0.25, 42);
    write_wav("guides/audio/white_noise.wav", buf, N);

    /* Impulse train — 4 clicks at 0.5s intervals */
    memset(buf, 0, N * sizeof(double));
    for (int i = 0; i < 4; i++) {
        unsigned pos = (unsigned)(i * 0.5 * SAMPLE_RATE);
        if (pos < N)
            buf[pos] = 0.8;
    }
    write_wav("guides/audio/impulse_train.wav", buf, N);

    /* --------------------------------------------------------------------
     * Simple effects: before/after clips
     * ------------------------------------------------------------------*/

    /* Delay source: percussive click train. */
    memset(buf, 0, N * sizeof(double));
    for (unsigned i = 0; i < N; i += (unsigned)(0.35 * SAMPLE_RATE)) {
        if (i < N) buf[i] = 0.9;
    }
    write_wav("guides/audio/effect_delay_before.wav", buf, N);
    MD_delay_echo(buf, fx, N, 11025, 0.45, 1.0, 0.6);
    write_wav("guides/audio/effect_delay_after.wav", fx, N);

    /* Tremolo source: steady 220 Hz sine. */
    MD_sine_wave(buf, N, 0.8, 220.0, SAMPLE_RATE);
    write_wav("guides/audio/effect_tremolo_before.wav", buf, N);
    MD_tremolo(buf, fx, N, 5.0, 0.8, SAMPLE_RATE);
    write_wav("guides/audio/effect_tremolo_after.wav", fx, N);

    /* Comb source: short decaying tone burst. */
    memset(buf, 0, N * sizeof(double));
    unsigned burst_len = SAMPLE_RATE / 5; /* 200 ms burst */
    for (unsigned i = 0; i < burst_len; i++) {
        double env = exp(-6.0 * (double)i / (double)burst_len);
        buf[i] = 0.85 * env * sin(2.0 * M_PI * 330.0 * (double)i / SAMPLE_RATE);
    }
    write_wav("guides/audio/effect_comb_before.wav", buf, N);
    MD_comb_reverb(buf, fx, N, 1323, 0.75, 0.7, 0.6);
    write_wav("guides/audio/effect_comb_after.wav", fx, N);

    free(fx);
    free(buf);
    return 0;
}
