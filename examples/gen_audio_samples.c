/**
 * @file gen_audio_samples.c
 * @brief Docs-only utility: generate WAV audio samples for the signal generators guide.
 *
 * This program generates 7 WAV files into guides/audio/ so that the Doxygen
 * documentation can embed playable audio previews of each signal generator.
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

static int write_wav(const char *path, const double *buf, unsigned N)
{
    float *fbuf = malloc(N * sizeof(float));
    if (!fbuf) {
        fprintf(stderr, "allocation failed for %s\n", path);
        return -1;
    }
    for (unsigned i = 0; i < N; i++)
        fbuf[i] = (float)buf[i];

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
    if (!buf) {
        fprintf(stderr, "allocation failed\n");
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

    free(buf);
    return 0;
}
