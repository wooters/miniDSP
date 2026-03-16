/**
 * @file test_steg.c
 * @brief Tests for MD_steg (audio steganography).
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * Tests for MD_steg (audio steganography)
 * -----------------------------------------------------------------------*/

/** LSB encode-then-decode round-trip recovers the original message. */
static int test_steg_lsb_roundtrip(void)
{
    const unsigned N = 8000;
    double *host   = malloc(N * sizeof(double));
    double *stego  = malloc(N * sizeof(double));

    /* Use a sine wave as the host signal. */
    MD_sine_wave(host, N, 0.8, 440.0, 44100.0);

    const char *secret = "Hello, miniDSP!";
    unsigned enc = MD_steg_encode(host, stego, N, 44100.0,
                                  secret, MD_STEG_LSB);
    if (enc == 0) { free(host); free(stego); return 0; }

    char recovered[256];
    unsigned dec = MD_steg_decode(stego, N, 44100.0,
                                  recovered, 256, MD_STEG_LSB);

    int ok = (dec == enc) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** LSB capacity matches expected formula: (signal_len - 32) / 8. */
static int test_steg_lsb_capacity(void)
{
    unsigned cap = MD_steg_capacity(8032, 44100.0, MD_STEG_LSB);
    /* (8032 - 32) / 8 = 1000 */
    return cap == 1000;
}

/** LSB encoding introduces negligible distortion (< -80 dB relative). */
static int test_steg_lsb_distortion(void)
{
    const unsigned N = 4000;
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.8, 440.0, 44100.0);
    MD_steg_encode(host, stego, N, 44100.0, "test", MD_STEG_LSB);

    /* Compute maximum absolute difference. */
    double max_diff = 0.0;
    for (unsigned i = 0; i < N; i++) {
        double d = fabs(host[i] - stego[i]);
        if (d > max_diff) max_diff = d;
    }

    free(stego);
    free(host);

    /* The double→PCM16→double roundtrip introduces one quantisation step,
     * and flipping the LSB adds at most one more.  Worst case ≈ 2/32767
     * ≈ 6.1e-5.  Use a comfortable margin. */
    return max_diff < 7.0e-5;
}

/** LSB handles longer messages that use more of the signal. */
static int test_steg_lsb_long_message(void)
{
    const unsigned N = 44100;  /* 1 second at 44.1 kHz */
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.5, 261.63, 44100.0);

    /* Build a long message. */
    char msg[512];
    memset(msg, 0, sizeof(msg));
    for (int i = 0; i < 500; i++)
        msg[i] = 'A' + (i % 26);

    unsigned enc = MD_steg_encode(host, stego, N, 44100.0,
                                  msg, MD_STEG_LSB);

    char recovered[512];
    unsigned dec = MD_steg_decode(stego, N, 44100.0,
                                  recovered, 512, MD_STEG_LSB);

    int ok = (dec == enc) && (dec > 0) && (strncmp(recovered, msg, dec) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Frequency-band encode-then-decode round-trip recovers the message. */
static int test_steg_freq_roundtrip(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 2.0);  /* 2 seconds */
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.5, 440.0, sr);

    const char *secret = "Hidden!";
    unsigned enc = MD_steg_encode(host, stego, N, sr,
                                  secret, MD_STEG_FREQ_BAND);
    if (enc == 0) { free(host); free(stego); return 0; }

    char recovered[256];
    unsigned dec = MD_steg_decode(stego, N, sr,
                                  recovered, 256, MD_STEG_FREQ_BAND);

    int ok = (dec == enc) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Frequency-band capacity is positive for a long enough signal at 44.1 kHz. */
static int test_steg_freq_capacity(void)
{
    unsigned cap = MD_steg_capacity(88200, 44100.0, MD_STEG_FREQ_BAND);
    /* chip_samples = 3.0 * 44100 / 1000 = 132.3 → 132
     * total_chips = 88200 / 132 = 668
     * (668 - 32) / 8 = 79 */
    return cap > 0 && cap <= 100;
}

/** Frequency-band decoding works with light additive noise. */
static int test_steg_freq_with_noise(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 2.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.5, 440.0, sr);

    const char *secret = "Noisy";
    unsigned enc = MD_steg_encode(host, stego, N, sr,
                                  secret, MD_STEG_FREQ_BAND);
    if (enc == 0) { free(host); free(stego); return 0; }

    /* Add low-level white noise (sigma=0.005). */
    double *noise = malloc(N * sizeof(double));
    MD_white_noise(noise, N, 0.005, 123);
    for (unsigned i = 0; i < N; i++)
        stego[i] += noise[i];
    free(noise);

    char recovered[256];
    unsigned dec = MD_steg_decode(stego, N, sr,
                                  recovered, 256, MD_STEG_FREQ_BAND);

    int ok = (dec == enc) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** LSB bytes roundtrip with embedded 0x00 bytes. */
static int test_steg_bytes_lsb_roundtrip(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.8, 440.0, sr);

    /* Data with embedded null bytes. */
    unsigned char data[] = {0x89, 0x50, 0x4E, 0x47, 0x00, 0x00, 0x0D, 0x0A,
                            0x1A, 0x0A, 0x00, 0x00, 0x00, 0x0D, 0xFF, 0xFE};
    unsigned data_len = sizeof(data);

    unsigned enc = MD_steg_encode_bytes(host, stego, N, sr,
                                        data, data_len, MD_STEG_LSB);
    if (enc != data_len) { free(host); free(stego); return 0; }

    unsigned char recovered[256];
    unsigned dec = MD_steg_decode_bytes(stego, N, sr,
                                        recovered, sizeof(recovered),
                                        MD_STEG_LSB);

    int ok = (dec == data_len) && (memcmp(recovered, data, data_len) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Frequency-band bytes roundtrip with embedded 0x00 bytes. */
static int test_steg_bytes_freq_roundtrip(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 2.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.5, 440.0, sr);

    unsigned char data[] = {0x00, 0x01, 0x02, 0x00, 0xFF, 0x00, 0xAB, 0xCD};
    unsigned data_len = sizeof(data);

    unsigned enc = MD_steg_encode_bytes(host, stego, N, sr,
                                        data, data_len, MD_STEG_FREQ_BAND);
    if (enc != data_len) { free(host); free(stego); return 0; }

    unsigned char recovered[256];
    unsigned dec = MD_steg_decode_bytes(stego, N, sr,
                                        recovered, sizeof(recovered),
                                        MD_STEG_FREQ_BAND);

    int ok = (dec == data_len) && (memcmp(recovered, data, data_len) == 0);
    free(stego);
    free(host);
    return ok;
}

/** String API and bytes API are interoperable. */
static int test_steg_bytes_string_compat(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));

    MD_sine_wave(host, N, 0.8, 440.0, sr);

    const char *secret = "Hello, steg!";
    unsigned slen = (unsigned)strlen(secret);

    /* Encode with string API, decode with bytes API. */
    unsigned enc = MD_steg_encode(host, stego, N, sr, secret, MD_STEG_LSB);
    if (enc != slen) { free(host); free(stego); return 0; }

    unsigned char bytes_out[256];
    unsigned dec = MD_steg_decode_bytes(stego, N, sr,
                                        bytes_out, sizeof(bytes_out),
                                        MD_STEG_LSB);
    if (dec != slen || memcmp(bytes_out, secret, slen) != 0) {
        free(host); free(stego); return 0;
    }

    /* Encode with bytes API, decode with string API. */
    unsigned enc2 = MD_steg_encode_bytes(host, stego, N, sr,
                                         (const unsigned char *)secret, slen,
                                         MD_STEG_LSB);
    if (enc2 != slen) { free(host); free(stego); return 0; }

    char str_out[256];
    unsigned dec2 = MD_steg_decode(stego, N, sr, str_out, sizeof(str_out),
                                   MD_STEG_LSB);
    int ok = (dec2 == slen) && (strcmp(str_out, secret) == 0);

    free(stego);
    free(host);
    return ok;
}

static int test_steg_detect_lsb(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    MD_steg_encode(host, stego, N, sr, "detect me", MD_STEG_LSB);

    int payload_type = -1;
    int method = MD_steg_detect(stego, N, sr, &payload_type);
    int ok = (method == MD_STEG_LSB) && (payload_type == MD_STEG_TYPE_TEXT);

    free(stego);
    free(host);
    return ok;
}

static int test_steg_detect_freq(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 3.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    MD_steg_encode(host, stego, N, sr, "bfsk msg", MD_STEG_FREQ_BAND);

    int payload_type = -1;
    int method = MD_steg_detect(stego, N, sr, &payload_type);
    int ok = (method == MD_STEG_FREQ_BAND) && (payload_type == MD_STEG_TYPE_TEXT);

    free(stego);
    free(host);
    return ok;
}

static int test_steg_detect_clean(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *signal = malloc(N * sizeof(double));
    MD_sine_wave(signal, N, 0.8, 440.0, sr);

    int method = MD_steg_detect(signal, N, sr, NULL);
    int ok = (method == -1);

    free(signal);
    return ok;
}

static int test_steg_detect_roundtrip(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    /* Encode binary data. */
    const unsigned char data[] = {0xDE, 0xAD, 0xBE, 0xEF, 0x00, 0x42};
    unsigned data_len = sizeof(data);
    MD_steg_encode_bytes(host, stego, N, sr, data, data_len, MD_STEG_LSB);

    /* Detect method and payload type. */
    int payload_type = -1;
    int method = MD_steg_detect(stego, N, sr, &payload_type);
    if (method != MD_STEG_LSB || payload_type != MD_STEG_TYPE_BINARY) {
        free(stego); free(host); return 0;
    }

    /* Decode with detected method and verify. */
    unsigned char recovered[256];
    unsigned decoded = MD_steg_decode_bytes(stego, N, sr,
                                            recovered, sizeof(recovered),
                                            method);
    int ok = (decoded == data_len) && (memcmp(recovered, data, data_len) == 0);

    free(stego);
    free(host);
    return ok;
}

/** Empty message should encode zero bits and leave host unchanged. */
static int test_steg_encode_empty_message(void)
{
    const unsigned N = 4000;
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, 44100.0);

    unsigned bits = MD_steg_encode(host, stego, N, 44100.0,
                                   "", MD_STEG_LSB);
    int ok = (bits == 0);

    free(stego);
    free(host);
    return ok;
}

/* -----------------------------------------------------------------------
 * Tests for MD_STEG_SPECTEXT (hybrid LSB + spectrogram text art)
 * -----------------------------------------------------------------------*/

/** Spectext capacity at 3 seconds / 44.1 kHz is limited by visual. */
static int test_steg_spectext_capacity(void)
{
    /* 3 sec at 44.1 kHz => visual cap = floor(3.0 / 0.24) = 12 */
    unsigned cap = MD_steg_capacity(132300, 44100.0, MD_STEG_SPECTEXT);
    return cap == 12;
}

/** Spectext capacity at 30 seconds / 48 kHz. */
static int test_steg_spectext_capacity_long(void)
{
    /* 30 sec => visual cap = floor(30.0 / 0.24) = 125 */
    unsigned cap = MD_steg_capacity(1440000, 48000.0, MD_STEG_SPECTEXT);
    return cap == 125;
}

/** Round-trip encode/decode at 48 kHz (no resampling needed). */
static int test_steg_spectext_roundtrip_48k(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 3.0);  /* 3 seconds */
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    const char *secret = "HELLO";
    unsigned enc = MD_steg_encode(host, stego, N, sr,
                                  secret, MD_STEG_SPECTEXT);
    if (enc == 0) { free(host); free(stego); return 0; }

    char recovered[256];
    unsigned dec = MD_steg_decode(stego, N, sr,
                                  recovered, 256, MD_STEG_SPECTEXT);

    int ok = (dec == enc) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Round-trip encode/decode at 44.1 kHz (resampling triggered). */
static int test_steg_spectext_roundtrip_44k(void)
{
    const double sr = 44100.0;
    const unsigned N = (unsigned)(sr * 3.0);  /* 3 seconds */
    double *host = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    /* Output is at 48 kHz — compute required length. */
    unsigned out_len = MD_resample_output_len(N, sr, 48000.0);
    double *stego = malloc(out_len * sizeof(double));

    const char *secret = "HI";
    unsigned enc = MD_steg_encode(host, stego, N, sr,
                                  secret, MD_STEG_SPECTEXT);
    if (enc == 0) { free(host); free(stego); return 0; }

    /* Decode from the 48 kHz output. */
    char recovered[256];
    unsigned dec = MD_steg_decode(stego, out_len, 48000.0,
                                  recovered, 256, MD_STEG_SPECTEXT);

    int ok = (dec == enc) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Visual truncation: 20-char message in 3-sec host, full message via decode. */
static int test_steg_spectext_truncation(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 3.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    /* 20 chars, but visual cap is 12 at 3 sec. LSB can hold all 20. */
    const char *secret = "TWENTY_CHAR_MESSAGE!";
    unsigned enc = MD_steg_encode(host, stego, N, sr,
                                  secret, MD_STEG_SPECTEXT);
    if (enc != 20) { free(host); free(stego); return 0; }

    char recovered[256];
    unsigned dec = MD_steg_decode(stego, N, sr,
                                  recovered, 256, MD_STEG_SPECTEXT);

    int ok = (dec == 20) && (strcmp(recovered, secret) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Binary encode/decode round-trip via spectext. */
static int test_steg_spectext_bytes_roundtrip(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 5.0);  /* 5 sec for capacity */
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.5, 440.0, sr);

    unsigned char data[] = {0x89, 0x50, 0x4E, 0x00, 0xFF, 0xAB};
    unsigned data_len = sizeof(data);

    unsigned enc = MD_steg_encode_bytes(host, stego, N, sr,
                                        data, data_len, MD_STEG_SPECTEXT);
    if (enc != data_len) { free(host); free(stego); return 0; }

    unsigned char recovered[256];
    unsigned dec = MD_steg_decode_bytes(stego, N, sr,
                                        recovered, sizeof(recovered),
                                        MD_STEG_SPECTEXT);

    int ok = (dec == data_len) && (memcmp(recovered, data, data_len) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Spectext decode produces same result as LSB decode on same signal. */
static int test_steg_spectext_decode_equals_lsb(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 3.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    MD_steg_encode(host, stego, N, sr, "MATCH", MD_STEG_SPECTEXT);

    char msg_spectext[256], msg_lsb[256];
    unsigned dec_st = MD_steg_decode(stego, N, sr,
                                     msg_spectext, 256, MD_STEG_SPECTEXT);
    unsigned dec_lsb = MD_steg_decode(stego, N, sr,
                                      msg_lsb, 256, MD_STEG_LSB);

    int ok = (dec_st == dec_lsb) && (strcmp(msg_spectext, msg_lsb) == 0);
    free(stego);
    free(host);
    return ok;
}

/** Detect spectext-encoded signal returns MD_STEG_SPECTEXT. */
static int test_steg_detect_spectext(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 3.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    MD_steg_encode(host, stego, N, sr, "VISIBLE", MD_STEG_SPECTEXT);

    int payload_type = -1;
    int method = MD_steg_detect(stego, N, sr, &payload_type);
    int ok = (method == MD_STEG_SPECTEXT) && (payload_type == MD_STEG_TYPE_TEXT);

    free(stego);
    free(host);
    return ok;
}

/** Detect plain LSB at 48 kHz still returns MD_STEG_LSB (no ultrasonic energy). */
static int test_steg_detect_lsb_not_spectext(void)
{
    const double sr = 48000.0;
    const unsigned N = (unsigned)(sr * 1.0);
    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    MD_sine_wave(host, N, 0.8, 440.0, sr);

    MD_steg_encode(host, stego, N, sr, "just lsb", MD_STEG_LSB);

    int payload_type = -1;
    int method = MD_steg_detect(stego, N, sr, &payload_type);
    int ok = (method == MD_STEG_LSB) && (payload_type == MD_STEG_TYPE_TEXT);

    free(stego);
    free(host);
    return ok;
}

/* -----------------------------------------------------------------------
 * Public entry point
 * -----------------------------------------------------------------------*/

void run_steg_tests(void)
{
    printf("\n--- MD_steg (audio steganography) ---\n");
    RUN_TEST(test_steg_lsb_roundtrip);
    RUN_TEST(test_steg_lsb_capacity);
    RUN_TEST(test_steg_lsb_distortion);
    RUN_TEST(test_steg_lsb_long_message);
    RUN_TEST(test_steg_freq_roundtrip);
    RUN_TEST(test_steg_freq_capacity);
    RUN_TEST(test_steg_freq_with_noise);
    RUN_TEST(test_steg_bytes_lsb_roundtrip);
    RUN_TEST(test_steg_bytes_freq_roundtrip);
    RUN_TEST(test_steg_bytes_string_compat);
    RUN_TEST(test_steg_detect_lsb);
    RUN_TEST(test_steg_detect_freq);
    RUN_TEST(test_steg_detect_clean);
    RUN_TEST(test_steg_detect_roundtrip);
    RUN_TEST(test_steg_encode_empty_message);
    RUN_TEST(test_steg_spectext_capacity);
    RUN_TEST(test_steg_spectext_capacity_long);
    RUN_TEST(test_steg_spectext_roundtrip_48k);
    RUN_TEST(test_steg_spectext_roundtrip_44k);
    RUN_TEST(test_steg_spectext_truncation);
    RUN_TEST(test_steg_spectext_bytes_roundtrip);
    RUN_TEST(test_steg_spectext_decode_equals_lsb);
    RUN_TEST(test_steg_detect_spectext);
    RUN_TEST(test_steg_detect_lsb_not_spectext);
}
