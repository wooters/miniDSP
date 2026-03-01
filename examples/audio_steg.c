/**
 * @file audio_steg.c
 * @brief Example: audio steganography — hide and recover secret messages.
 *
 * Three modes of operation:
 *
 *   1) No arguments — self-test: generate a host signal, encode a message
 *      using both LSB and frequency-band methods, decode, and verify.
 *
 *   2) --encode METHOD MESSAGE [-i HOST.wav] [-o STEGO.wav]
 *      Encode a message into a host WAV file.  METHOD is "lsb" or "freq".
 *      If no host is given, a 3-second 440 Hz sine at 44.1 kHz is used.
 *
 *   3) --decode METHOD FILE
 *      Decode a hidden message from a stego WAV file.
 *
 * Build and run (from the repository root):
 *   make -C examples audio_steg
 *   cd examples && ./audio_steg
 *   cd examples && ./audio_steg --encode lsb "secret message" -o stego.wav
 *   cd examples && ./audio_steg --decode lsb stego.wav
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minidsp.h"
#include "fileio.h"

/* --------------------------------------------------------------------- */

static int parse_method(const char *str)
{
    if (strcmp(str, "lsb") == 0)  return MD_STEG_LSB;
    if (strcmp(str, "freq") == 0) return MD_STEG_FREQ_BAND;
    return -1;
}

static const char *method_name(int method)
{
    return (method == MD_STEG_LSB) ? "LSB" : "frequency-band";
}

/* --------------------------------------------------------------------- */

//! [encode-wav]
static int encode_wav(int method, const char *message,
                      const char *infile, const char *outfile)
{
    double  *host    = nullptr;
    unsigned signal_len;
    unsigned samprate;

    if (infile) {
        /* Read host from WAV. */
        float  *fdata   = nullptr;
        size_t  datalen = 0;
        unsigned sr     = 0;

        if (FIO_read_audio(infile, &fdata, &datalen, &sr, 1) != 0) {
            fprintf(stderr, "Error reading %s\n", infile);
            return 1;
        }
        if (datalen == 0 || datalen > UINT_MAX) {
            fprintf(stderr, "Invalid audio data in %s\n", infile);
            free(fdata);
            return 1;
        }
        signal_len = (unsigned)datalen;
        samprate   = sr;
        host = malloc(signal_len * sizeof(double));
        if (!host) { free(fdata); return 1; }
        for (unsigned i = 0; i < signal_len; i++)
            host[i] = (double)fdata[i];
        free(fdata);
    } else {
        /* Generate a default host signal: 3 s sine at 44.1 kHz. */
        samprate   = 44100;
        signal_len = samprate * 3;
        host = malloc(signal_len * sizeof(double));
        if (!host) return 1;
        MD_sine_wave(host, signal_len, 0.8, 440.0, (double)samprate);
        printf("No host file specified; using 3 s, 440 Hz sine at %u Hz\n",
               samprate);
    }

    if (method == MD_STEG_FREQ_BAND && samprate < 40000) {
        fprintf(stderr,
            "Frequency-band method requires sample rate >= 40 kHz "
            "(file is %u Hz).\n"
            "Use a 44.1 kHz or 48 kHz host, or use the LSB method.\n",
            samprate);
        free(host);
        return 1;
    }

    unsigned capacity = MD_steg_capacity(signal_len, (double)samprate, method);
    unsigned msg_len  = (unsigned)strlen(message);
    printf("Method: %s | Capacity: %u bytes | Message: %u bytes\n",
           method_name(method), capacity, msg_len);

    if (msg_len > capacity) {
        fprintf(stderr,
            "Message too long (%u bytes) for host signal capacity (%u bytes).\n"
            "Use a longer host signal or a shorter message.\n",
            msg_len, capacity);
        free(host);
        return 1;
    }

    double *stego = malloc(signal_len * sizeof(double));
    if (!stego) { free(host); return 1; }

    unsigned encoded = MD_steg_encode(host, stego, signal_len,
                                      (double)samprate, message, method);

    /* Write to WAV. */
    float *fout = malloc(signal_len * sizeof(float));
    if (!fout) { free(stego); free(host); return 1; }
    for (unsigned i = 0; i < signal_len; i++)
        fout[i] = (float)stego[i];

    int ret = FIO_write_wav(outfile, fout, signal_len, samprate);
    if (ret == 0)
        printf("Encoded %u bytes -> %s  (%u samples, %.3f s)\n",
               encoded, outfile, signal_len,
               (double)signal_len / (double)samprate);
    else
        fprintf(stderr, "Error writing %s\n", outfile);

    free(fout);
    free(stego);
    free(host);
    return ret;
}
//! [encode-wav]

/* --------------------------------------------------------------------- */

//! [decode-wav]
static int decode_wav(int method, const char *infile)
{
    float  *fdata   = nullptr;
    size_t  datalen = 0;
    unsigned samprate = 0;

    if (FIO_read_audio(infile, &fdata, &datalen, &samprate, 1) != 0) {
        fprintf(stderr, "Error reading %s\n", infile);
        return 1;
    }
    if (datalen == 0 || datalen > UINT_MAX) {
        fprintf(stderr, "Invalid audio in %s\n", infile);
        free(fdata);
        return 1;
    }

    printf("Read %s: %zu samples at %u Hz (%.3f s)\n",
           infile, datalen, samprate,
           (double)datalen / (double)samprate);

    double *stego = malloc(datalen * sizeof(double));
    if (!stego) { free(fdata); return 1; }
    for (size_t i = 0; i < datalen; i++)
        stego[i] = (double)fdata[i];
    free(fdata);

    char message[4096];
    unsigned decoded = MD_steg_decode(stego, (unsigned)datalen,
                                      (double)samprate,
                                      message, sizeof(message), method);
    free(stego);

    if (decoded == 0) {
        printf("No hidden message found (method: %s).\n", method_name(method));
        return 1;
    }

    printf("\nDecoded %u bytes (method: %s):\n  \"%s\"\n",
           decoded, method_name(method), message);
    return 0;
}
//! [decode-wav]

/* --------------------------------------------------------------------- */

//! [self-test]
static int self_test(void)
{
    const double   sr = 44100.0;
    const unsigned N  = (unsigned)(sr * 3.0);  /* 3 seconds */
    const char    *secret = "The quick brown fox jumps over the lazy dog.";

    double *host  = malloc(N * sizeof(double));
    double *stego = malloc(N * sizeof(double));
    char    recovered[256];

    if (!host || !stego) {
        fprintf(stderr, "allocation failed\n");
        free(stego); free(host);
        return 1;
    }

    MD_sine_wave(host, N, 0.8, 440.0, sr);

    int pass = 1;

    /* --- LSB test --- */
    printf("=== LSB steganography ===\n");
    printf("  Host: 3 s sine wave at 440 Hz, %.0f Hz sample rate\n", sr);
    printf("  Capacity: %u bytes\n",
           MD_steg_capacity(N, sr, MD_STEG_LSB));
    printf("  Message (%zu bytes): \"%s\"\n", strlen(secret), secret);

    unsigned enc_lsb = MD_steg_encode(host, stego, N, sr,
                                      secret, MD_STEG_LSB);
    printf("  Encoded: %u bytes\n", enc_lsb);

    unsigned dec_lsb = MD_steg_decode(stego, N, sr,
                                      recovered, sizeof(recovered),
                                      MD_STEG_LSB);
    printf("  Decoded: %u bytes -> \"%s\"\n", dec_lsb, recovered);

    if (dec_lsb != enc_lsb || strcmp(recovered, secret) != 0) {
        printf("  LSB FAILED: decoded message does not match!\n");
        pass = 0;
    } else {
        /* Compute distortion. */
        double max_diff = 0.0;
        for (unsigned i = 0; i < N; i++) {
            double d = fabs(host[i] - stego[i]);
            if (d > max_diff) max_diff = d;
        }
        printf("  Max distortion: %.2e (%.1f dB)\n",
               max_diff, 20.0 * log10(max_diff + 1e-30));
        printf("  LSB PASSED\n");
    }

    /* --- Frequency-band test --- */
    printf("\n=== Frequency-band steganography (BFSK) ===\n");
    printf("  Host: 3 s sine wave at 440 Hz, %.0f Hz sample rate\n", sr);
    printf("  Capacity: %u bytes\n",
           MD_steg_capacity(N, sr, MD_STEG_FREQ_BAND));
    printf("  Message (%zu bytes): \"%s\"\n", strlen(secret), secret);

    unsigned enc_freq = MD_steg_encode(host, stego, N, sr,
                                       secret, MD_STEG_FREQ_BAND);
    printf("  Encoded: %u bytes\n", enc_freq);

    unsigned dec_freq = MD_steg_decode(stego, N, sr,
                                       recovered, sizeof(recovered),
                                       MD_STEG_FREQ_BAND);
    printf("  Decoded: %u bytes -> \"%s\"\n", dec_freq, recovered);

    if (dec_freq != enc_freq || strcmp(recovered, secret) != 0) {
        printf("  Frequency-band FAILED: decoded message does not match!\n");
        pass = 0;
    } else {
        printf("  Frequency-band PASSED\n");
    }

    if (pass)
        printf("\nSelf-test PASSED: both methods recovered the message.\n");
    else
        printf("\nSelf-test FAILED.\n");

    free(stego);
    free(host);
    return pass ? 0 : 1;
}
//! [self-test]

/* --------------------------------------------------------------------- */

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  %s                                    (self-test)\n"
        "  %s --encode METHOD MSG [-i HOST] [-o OUT]\n"
        "  %s --decode METHOD FILE\n"
        "\n"
        "METHOD: \"lsb\" or \"freq\"\n"
        "\n"
        "Examples:\n"
        "  %s --encode lsb \"secret message\" -o stego.wav\n"
        "  %s --decode lsb stego.wav\n"
        "  %s --encode freq \"hidden\" -i music.wav -o stego.wav\n"
        "  %s --decode freq stego.wav\n",
        prog, prog, prog, prog, prog, prog, prog);
}

int main(int argc, char *argv[])
{
    if (argc == 1)
        return self_test();

    if (argc >= 4 && strcmp(argv[1], "--encode") == 0) {
        int method = parse_method(argv[2]);
        if (method < 0) {
            fprintf(stderr, "Unknown method '%s' (use 'lsb' or 'freq')\n",
                    argv[2]);
            return 1;
        }
        const char *message = argv[3];
        const char *infile  = nullptr;
        const char *outfile = "steg_output.wav";

        for (int i = 4; i < argc - 1; i++) {
            if (strcmp(argv[i], "-i") == 0) infile  = argv[++i];
            else if (strcmp(argv[i], "-o") == 0) outfile = argv[++i];
        }
        return encode_wav(method, message, infile, outfile);
    }

    if (argc >= 4 && strcmp(argv[1], "--decode") == 0) {
        int method = parse_method(argv[2]);
        if (method < 0) {
            fprintf(stderr, "Unknown method '%s' (use 'lsb' or 'freq')\n",
                    argv[2]);
            return 1;
        }
        return decode_wav(method, argv[3]);
    }

    usage(argv[0]);
    return 1;
}
