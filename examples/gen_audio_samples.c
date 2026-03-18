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

static int write_wav_sr(const char *path, const double *buf, unsigned N,
                         unsigned sr)
{
    float *fbuf = malloc(N * sizeof(float));
    if (!fbuf) {
        fprintf(stderr, "allocation failed for %s\n", path);
        return -1;
    }
    for (unsigned i = 0; i < N; i++)
        fbuf[i] = (float)clamp_unit(buf[i]);

    int rc = FIO_write_wav(path, fbuf, N, sr);
    free(fbuf);
    if (rc != 0)
        fprintf(stderr, "failed to write %s\n", path);
    else
        printf("  %s\n", path);
    return rc;
}

int main(void)
{
    MD_set_error_handler(NULL);  /* use default stderr handler */

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

    /* --------------------------------------------------------------------
     * Spectrogram text: "HELLO" at 16 kHz
     * ------------------------------------------------------------------*/
    {
        const unsigned st_sr = 16000;
        const double st_dur = 2.25;
        const double st_pad = 0.5;   /* silence before and after text */
        unsigned st_max = (unsigned)(st_sr * st_dur) + 1024;
        double *st_text = malloc(st_max * sizeof(double));
        if (!st_text) {
            fprintf(stderr, "allocation failed for spectrogram text\n");
        } else {
            unsigned st_text_n = MD_spectrogram_text(st_text, st_max, "HELLO",
                                                     400.0, 7300.0, st_dur, st_sr);
            /* Pad with silence before and after */
            unsigned st_pad_samp = (unsigned)(st_pad * st_sr);
            unsigned st_n = st_pad_samp + st_text_n + st_pad_samp;
            double *st_padded = calloc(st_n, sizeof(double));
            if (!st_padded) {
                fprintf(stderr, "allocation failed for padded spectrogram text\n");
            } else {
                memcpy(st_padded + st_pad_samp, st_text, st_text_n * sizeof(double));
                /* Write at 16 kHz — inline FIO_write_wav since write_wav() hardcodes 44100 */
                float *st_fbuf = malloc(st_n * sizeof(float));
                if (!st_fbuf) {
                    fprintf(stderr, "allocation failed for spectrogram text float buf\n");
                } else {
                    for (unsigned i = 0; i < st_n; i++)
                        st_fbuf[i] = (float)clamp_unit(st_padded[i]);
                    int rc = FIO_write_wav("guides/audio/spectrogram_text_hello.wav",
                                           st_fbuf, st_n, st_sr);
                    if (rc != 0)
                        fprintf(stderr, "failed to write spectrogram_text_hello.wav\n");
                    else
                        printf("  guides/audio/spectrogram_text_hello.wav\n");
                    free(st_fbuf);
                }
                free(st_padded);
            }
            free(st_text);
        }
    }

    /* --------------------------------------------------------------------
     * Shepard tones: rising, falling, and static (5 seconds at 44100 Hz)
     * ------------------------------------------------------------------*/
    {
        const unsigned shep_n = (unsigned)(SAMPLE_RATE * 5.0);
        double *shep = malloc(shep_n * sizeof(double));
        if (!shep) {
            fprintf(stderr, "allocation failed for shepard tone\n");
        } else {
            /* Rising */
            MD_shepard_tone(shep, shep_n, 0.8, 440.0, SAMPLE_RATE, 0.5, 8);
            write_wav("guides/audio/shepard_rising.wav", shep, shep_n);

            /* Falling */
            MD_shepard_tone(shep, shep_n, 0.8, 440.0, SAMPLE_RATE, -0.5, 8);
            write_wav("guides/audio/shepard_falling.wav", shep, shep_n);

            /* Static chord */
            MD_shepard_tone(shep, shep_n, 0.8, 440.0, SAMPLE_RATE, 0.0, 8);
            write_wav("guides/audio/shepard_static.wav", shep, shep_n);

            free(shep);
        }
    }

    /* --------------------------------------------------------------------
     * Audio steganography: host, LSB stego, and freq-band stego
     * ------------------------------------------------------------------*/
    {
        const unsigned steg_n = (unsigned)(SAMPLE_RATE * 3.0);
        double *steg_host  = malloc(steg_n * sizeof(double));
        double *steg_out   = malloc(steg_n * sizeof(double));
        if (!steg_host || !steg_out) {
            fprintf(stderr, "allocation failed for steganography samples\n");
        } else {
            const char *secret = "Hidden message inside audio!";

            /* Host: 440 Hz sine, 3 seconds at 44100 Hz */
            MD_sine_wave(steg_host, steg_n, 0.8, 440.0, SAMPLE_RATE);
            write_wav("guides/audio/steg_host.wav", steg_host, steg_n);

            /* LSB stego */
            MD_steg_encode(steg_host, steg_out, steg_n, SAMPLE_RATE,
                           secret, MD_STEG_LSB);
            write_wav("guides/audio/steg_lsb.wav", steg_out, steg_n);

            /* Frequency-band stego */
            MD_steg_encode(steg_host, steg_out, steg_n, SAMPLE_RATE,
                           secret, MD_STEG_FREQ_BAND);
            write_wav("guides/audio/steg_freq.wav", steg_out, steg_n);

            free(steg_out);
            free(steg_host);
        }
    }

    /* --------------------------------------------------------------------
     * Spectrogram-text steganography: sa2.wav host with "miniDSP" at 48 kHz
     * ------------------------------------------------------------------*/
    {
        float *fdata = NULL;
        size_t datalen = 0;
        unsigned sr = 0;
        if (FIO_read_audio("samples/sa2.wav", &fdata, &datalen, &sr, 1) != 0) {
            fprintf(stderr, "failed to read samples/sa2.wav for spectext\n");
        } else {
            unsigned host_n = (unsigned)datalen;
            double *host = malloc(host_n * sizeof(double));
            if (!host) {
                fprintf(stderr, "allocation failed for spectext host\n");
            } else {
                for (unsigned i = 0; i < host_n; i++)
                    host[i] = (double)fdata[i];

                unsigned out_n = MD_resample_output_len(host_n, (double)sr,
                                                         48000.0);
                double *steg_spec = malloc(out_n * sizeof(double));
                if (!steg_spec) {
                    fprintf(stderr, "allocation failed for spectext output\n");
                } else {
                    MD_steg_encode(host, steg_spec, host_n, (double)sr,
                                   "miniDSP", MD_STEG_SPECTEXT);
                    write_wav_sr("guides/audio/steg_spectext.wav",
                                 steg_spec, out_n, 48000);
                    free(steg_spec);
                }
                free(host);
            }
            free(fdata);
        }
    }

    /* --------------------------------------------------------------------
     * Binary steganography: hide PNG images in audio, recover them
     * ------------------------------------------------------------------*/
    {
        const char *images[] = {
            "tools/audio_steg/space_invader.png", "tools/audio_steg/minidsp_qr.png"
        };
        const char *wav_out[] = {
            "guides/audio/steg_lsb_invader.wav",
            "guides/audio/steg_lsb_qr.wav"
        };
        const char *png_out[] = {
            "guides/images/steg_recovered_invader.png",
            "guides/images/steg_recovered_qr.png"
        };

        for (int img_i = 0; img_i < 2; img_i++) {
            /* Read the source PNG */
            FILE *fp = fopen(images[img_i], "rb");
            if (!fp) {
                fprintf(stderr, "cannot open %s\n", images[img_i]);
                continue;
            }
            fseek(fp, 0, SEEK_END);
            unsigned data_len = (unsigned)ftell(fp);
            fseek(fp, 0, SEEK_SET);
            unsigned char *img_data = malloc(data_len);
            if (!img_data) {
                fprintf(stderr, "allocation failed for %s\n", images[img_i]);
                fclose(fp);
                continue;
            }
            fread(img_data, 1, data_len, fp);
            fclose(fp);

            /* Generate host signal sized for the payload */
            unsigned sig_n = data_len * 8 + 32 + 1024;  /* payload + header + margin */
            double *host  = malloc(sig_n * sizeof(double));
            double *stego = malloc(sig_n * sizeof(double));
            if (!host || !stego) {
                fprintf(stderr, "allocation failed for steg buffers\n");
                free(stego);
                free(host);
                free(img_data);
                continue;
            }

            MD_sine_wave(host, sig_n, 0.8, 440.0, SAMPLE_RATE);

            /* Encode image into audio */
            unsigned encoded = MD_steg_encode_bytes(host, stego, sig_n,
                                                    (double)SAMPLE_RATE,
                                                    img_data, data_len,
                                                    MD_STEG_LSB);
            if (encoded == 0) {
                fprintf(stderr, "encode failed for %s\n", images[img_i]);
                free(stego);
                free(host);
                free(img_data);
                continue;
            }

            /* Write stego WAV */
            write_wav(wav_out[img_i], stego, sig_n);

            /* Decode to recover the image bytes */
            unsigned char *recovered = malloc(data_len);
            if (!recovered) {
                fprintf(stderr, "allocation failed for recovery buffer\n");
            } else {
                unsigned decoded = MD_steg_decode_bytes(stego, sig_n,
                                                        (double)SAMPLE_RATE,
                                                        recovered, data_len,
                                                        MD_STEG_LSB);
                if (decoded != data_len) {
                    fprintf(stderr, "decode length mismatch for %s: %u vs %u\n",
                            images[img_i], decoded, data_len);
                } else {
                    /* Write recovered PNG */
                    FILE *out = fopen(png_out[img_i], "wb");
                    if (!out) {
                        fprintf(stderr, "cannot write %s\n", png_out[img_i]);
                    } else {
                        fwrite(recovered, 1, data_len, out);
                        fclose(out);
                        printf("  %s\n", png_out[img_i]);
                    }
                }
                free(recovered);
            }

            free(stego);
            free(host);
            free(img_data);
        }
    }

    free(fx);
    free(buf);
    return 0;
}
