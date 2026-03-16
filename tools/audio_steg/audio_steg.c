/**
 * @file audio_steg.c
 * @brief Example: audio steganography — hide and recover secret messages or binary data.
 *
 * Six modes of operation:
 *
 *   1) No arguments — self-test: generate a host signal, encode a message
 *      using both LSB and frequency-band methods, decode, and verify.
 *
 *   2) --encode METHOD MESSAGE [-i HOST.wav] [-o STEGO.wav]
 *      Encode a text message into a host WAV file.
 *
 *   3) --decode [METHOD] FILE
 *      Decode a hidden text message from a stego WAV file.
 *      METHOD is optional — omit it to auto-detect.
 *
 *   4) --encode-image METHOD IMAGE [-i HOST.wav] [-o STEGO.wav]
 *      Encode a binary file (e.g. PNG image) into a host WAV file.
 *
 *   5) --decode-image [METHOD] FILE -o IMAGE_OUT
 *      Decode a binary file from a stego WAV file.
 *      METHOD is optional — omit it to auto-detect.
 *
 *   6) FILE (bare filename, no flags)
 *      Auto-detect method and payload type.  Text payloads are printed;
 *      binary payloads report size and suggest --decode-image.
 *
 * Build and run (from the repository root):
 *   make -C examples audio_steg
 *   cd examples && ./audio_steg
 *   cd examples && ./audio_steg --encode lsb "secret message" -o stego.wav
 *   cd examples && ./audio_steg --decode stego.wav
 *   cd examples && ./audio_steg stego.wav
 *   cd examples && ./audio_steg --encode-image lsb space_invader.png -o stego.wav
 *   cd examples && ./audio_steg --decode-image lsb stego.wav -o recovered.png
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minidsp.h"
#include "fileio.h"

/** Bits for the message-length header (matches library constant). */
#define HEADER_BITS 32

/* --------------------------------------------------------------------- */

/** Read a WAV file and convert to double.  Caller must free *out. */
static int read_wav_to_double(const char *path, double **out,
                              unsigned *out_len, unsigned *out_sr)
{
    float  *fdata   = NULL;
    size_t  datalen = 0;
    unsigned sr     = 0;

    if (FIO_read_audio(path, &fdata, &datalen, &sr, 1) != 0) {
        fprintf(stderr, "Error reading %s\n", path);
        return 1;
    }
    if (datalen == 0 || datalen > UINT_MAX) {
        fprintf(stderr, "Invalid audio data in %s\n", path);
        free(fdata);
        return 1;
    }

    *out_len = (unsigned)datalen;
    *out_sr  = sr;
    *out = malloc(datalen * sizeof(double));
    if (!*out) { free(fdata); return 1; }
    for (size_t i = 0; i < datalen; i++)
        (*out)[i] = (double)fdata[i];
    free(fdata);
    return 0;
}

/** Convert double signal to float and write as WAV. */
static int write_double_as_wav(const char *path, const double *signal,
                               unsigned len, unsigned samprate)
{
    float *fdata = malloc(len * sizeof(float));
    if (!fdata) return 1;
    for (unsigned i = 0; i < len; i++)
        fdata[i] = (float)signal[i];
    int ret = FIO_write_wav(path, fdata, len, samprate);
    free(fdata);
    return ret;
}

/* --------------------------------------------------------------------- */

static int parse_method(const char *str)
{
    if (strcmp(str, "lsb") == 0)      return MD_STEG_LSB;
    if (strcmp(str, "freq") == 0)     return MD_STEG_FREQ_BAND;
    if (strcmp(str, "spectext") == 0) return MD_STEG_SPECTEXT;
    return -1;
}

static const char *method_name(int method)
{
    if (method == MD_STEG_LSB)       return "LSB";
    if (method == MD_STEG_SPECTEXT)  return "spectext";
    return "frequency-band";
}

/** Return the CLI method string for --decode-image suggestions. */
static const char *method_cli_name(int method)
{
    if (method == MD_STEG_LSB)       return "lsb";
    if (method == MD_STEG_SPECTEXT)  return "spectext";
    return "freq";
}

/* --------------------------------------------------------------------- */

//! [encode-wav]
static int encode_wav(int method, const char *message,
                      const char *infile, const char *outfile)
{
    double  *host    = NULL;
    unsigned signal_len;
    unsigned samprate;

    if (infile) {
        if (read_wav_to_double(infile, &host, &signal_len, &samprate) != 0)
            return 1;
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

    /* Spectext outputs at 48 kHz — allocate for the larger signal. */
    unsigned out_len = signal_len;
    unsigned out_sr  = samprate;
    if (method == MD_STEG_SPECTEXT) {
        if (samprate < 48000) {
            out_len = MD_resample_output_len(signal_len, (double)samprate,
                                              48000.0);
        }
        out_sr = 48000;
    }

    double *stego = malloc(out_len * sizeof(double));
    if (!stego) { free(host); return 1; }

    unsigned encoded = MD_steg_encode(host, stego, signal_len,
                                      (double)samprate, message, method);

    int ret = write_double_as_wav(outfile, stego, out_len, out_sr);
    if (ret == 0)
        printf("Encoded %u bytes -> %s  (%u samples, %.3f s, %u Hz)\n",
               encoded, outfile, out_len,
               (double)out_len / (double)out_sr, out_sr);
    else
        fprintf(stderr, "Error writing %s\n", outfile);

    free(stego);
    free(host);
    return ret;
}
//! [encode-wav]

/* --------------------------------------------------------------------- */

//! [decode-wav]
static int decode_wav(int method, const char *infile)
{
    double  *stego    = NULL;
    unsigned signal_len;
    unsigned samprate;

    if (read_wav_to_double(infile, &stego, &signal_len, &samprate) != 0)
        return 1;

    printf("Read %s: %u samples at %u Hz (%.3f s)\n",
           infile, signal_len, samprate,
           (double)signal_len / (double)samprate);

    char message[4096];
    unsigned decoded = MD_steg_decode(stego, signal_len, (double)samprate,
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

//! [encode-image]
static int encode_image_wav(int method, const char *image_path,
                            const char *infile, const char *outfile)
{
    /* Read the image file into memory. */
    FILE *fp = fopen(image_path, "rb");
    if (!fp) {
        fprintf(stderr, "Cannot open image file: %s\n", image_path);
        return 1;
    }
    fseek(fp, 0, SEEK_END);
    long fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (fsize <= 0 || (unsigned long)fsize > UINT_MAX) {
        fprintf(stderr, "Invalid file size: %s\n", image_path);
        fclose(fp);
        return 1;
    }
    unsigned data_len = (unsigned)fsize;
    unsigned char *data = malloc(data_len);
    if (!data) { fclose(fp); return 1; }
    if (fread(data, 1, data_len, fp) != data_len) {
        fprintf(stderr, "Failed to read %s\n", image_path);
        free(data);
        fclose(fp);
        return 1;
    }
    fclose(fp);

    printf("Image: %s (%u bytes)\n", image_path, data_len);

    double  *host    = NULL;
    unsigned signal_len;
    unsigned samprate;

    if (infile) {
        if (read_wav_to_double(infile, &host, &signal_len, &samprate) != 0) {
            free(data);
            return 1;
        }
    } else {
        /* Generate a host signal sized for the payload. */
        samprate = 44100;
        unsigned min_samples = data_len * 8 + HEADER_BITS;
        if (method == MD_STEG_FREQ_BAND) {
            unsigned cs = (unsigned)(3.0 * samprate / 1000.0);
            min_samples = (data_len * 8 + HEADER_BITS) * cs;
        }
        unsigned half_sec = samprate / 2;
        if (min_samples < half_sec)
            min_samples = half_sec;
        signal_len = min_samples + min_samples / 10;  /* 10% margin */
        host = malloc(signal_len * sizeof(double));
        if (!host) { free(data); return 1; }
        MD_sine_wave(host, signal_len, 0.8, 440.0, (double)samprate);
        printf("Generated host: %u samples (%.3f s) at %u Hz\n",
               signal_len, (double)signal_len / (double)samprate, samprate);
    }

    if (method == MD_STEG_FREQ_BAND && samprate < 40000) {
        fprintf(stderr,
            "Frequency-band method requires sample rate >= 40 kHz "
            "(got %u Hz).\n", samprate);
        free(host);
        free(data);
        return 1;
    }

    unsigned capacity = MD_steg_capacity(signal_len, (double)samprate, method);
    printf("Method: %s | Capacity: %u bytes | Payload: %u bytes\n",
           method_name(method), capacity, data_len);

    if (data_len > capacity) {
        fprintf(stderr,
            "Payload too large (%u bytes) for host capacity (%u bytes).\n",
            data_len, capacity);
        free(host);
        free(data);
        return 1;
    }

    /* Spectext outputs at 48 kHz — allocate for the larger signal. */
    unsigned out_len = signal_len;
    unsigned out_sr  = samprate;
    if (method == MD_STEG_SPECTEXT) {
        if (samprate < 48000) {
            out_len = MD_resample_output_len(signal_len, (double)samprate,
                                              48000.0);
        }
        out_sr = 48000;
    }

    double *stego = malloc(out_len * sizeof(double));
    if (!stego) { free(host); free(data); return 1; }

    unsigned encoded = MD_steg_encode_bytes(host, stego, signal_len,
                                            (double)samprate,
                                            data, data_len, method);

    int ret = write_double_as_wav(outfile, stego, out_len, out_sr);
    if (ret == 0)
        printf("Encoded %u bytes -> %s  (%u samples, %.3f s, %u Hz)\n",
               encoded, outfile, out_len,
               (double)out_len / (double)out_sr, out_sr);
    else
        fprintf(stderr, "Error writing %s\n", outfile);

    free(stego);
    free(host);
    free(data);
    return ret;
}
//! [encode-image]

/* --------------------------------------------------------------------- */

//! [decode-image]
static int decode_image_wav(int method, const char *infile,
                            const char *outfile)
{
    double  *stego    = NULL;
    unsigned signal_len;
    unsigned samprate;

    if (read_wav_to_double(infile, &stego, &signal_len, &samprate) != 0)
        return 1;

    printf("Read %s: %u samples at %u Hz (%.3f s)\n",
           infile, signal_len, samprate,
           (double)signal_len / (double)samprate);

    /* Allocate a buffer large enough for the maximum capacity (capped at
     * 16 MB to avoid unbounded allocation from corrupted headers). */
    unsigned capacity = MD_steg_capacity(signal_len, (double)samprate, method);
    if (capacity > 16u * 1024 * 1024)
        capacity = 16u * 1024 * 1024;
    unsigned char *buf = malloc(capacity > 0 ? capacity : 1);
    if (!buf) { free(stego); return 1; }

    unsigned decoded = MD_steg_decode_bytes(stego, signal_len,
                                            (double)samprate,
                                            buf, capacity, method);
    free(stego);

    if (decoded == 0) {
        printf("No hidden data found (method: %s).\n", method_name(method));
        free(buf);
        return 1;
    }

    printf("Decoded %u bytes (method: %s)\n", decoded, method_name(method));

    FILE *fp = fopen(outfile, "wb");
    if (!fp) {
        fprintf(stderr, "Cannot open output file: %s\n", outfile);
        free(buf);
        return 1;
    }
    fwrite(buf, 1, decoded, fp);
    fclose(fp);
    free(buf);

    printf("Written to %s\n", outfile);
    return 0;
}
//! [decode-image]

/* --------------------------------------------------------------------- */

//! [auto-detect]
static int auto_decode_wav(const char *infile)
{
    double  *stego    = NULL;
    unsigned signal_len;
    unsigned samprate;

    if (read_wav_to_double(infile, &stego, &signal_len, &samprate) != 0)
        return 1;

    printf("Read %s: %u samples at %u Hz (%.3f s)\n",
           infile, signal_len, samprate,
           (double)signal_len / (double)samprate);

    int payload_type = -1;
    int method = MD_steg_detect(stego, signal_len, (double)samprate,
                                &payload_type);

    if (method < 0) {
        printf("No hidden payload detected.\n");
        free(stego);
        return 1;
    }

    printf("Detected: %s method, %s payload\n",
           method_name(method),
           payload_type == MD_STEG_TYPE_BINARY ? "binary" : "text");

    if (payload_type == MD_STEG_TYPE_BINARY) {
        /* Probe the payload size without a full decode. */
        unsigned capacity = MD_steg_capacity(signal_len, (double)samprate,
                                             method);
        if (capacity > 16u * 1024 * 1024)
            capacity = 16u * 1024 * 1024;
        unsigned char *buf = malloc(capacity > 0 ? capacity : 1);
        if (!buf) { free(stego); return 1; }

        unsigned decoded = MD_steg_decode_bytes(stego, signal_len,
                                                (double)samprate,
                                                buf, capacity, method);
        free(buf);
        free(stego);

        printf("Binary payload: %u bytes\n", decoded);
        printf("Use --decode-image %s %s -o OUTPUT to extract.\n",
               method_cli_name(method), infile);
        return 0;
    }

    /* Text payload — decode and print. */
    char message[4096];
    unsigned decoded = MD_steg_decode(stego, signal_len, (double)samprate,
                                      message, sizeof(message), method);
    free(stego);

    if (decoded == 0) {
        printf("Decode returned 0 bytes.\n");
        return 1;
    }

    printf("\nDecoded %u bytes:\n  \"%s\"\n", decoded, message);
    return 0;
}
//! [auto-detect]

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

    /* --- Spectext test (uses 48 kHz output) --- */
    printf("\n=== Spectrogram text steganography (spectext) ===\n");
    const char *spec_secret = "miniDSP";
    printf("  Host: 3 s sine wave at 440 Hz, %.0f Hz sample rate\n", sr);
    printf("  Capacity: %u chars\n",
           MD_steg_capacity(N, sr, MD_STEG_SPECTEXT));
    printf("  Message (%zu bytes): \"%s\"\n", strlen(spec_secret), spec_secret);

    /* Spectext may upsample to 48 kHz — allocate for larger output. */
    unsigned spec_out_len = MD_resample_output_len(N, sr, 48000.0);
    double *stego_spec = malloc(spec_out_len * sizeof(double));
    if (!stego_spec) {
        fprintf(stderr, "allocation failed\n");
        free(stego); free(host);
        return 1;
    }

    unsigned enc_spec = MD_steg_encode(host, stego_spec, N, sr,
                                        spec_secret, MD_STEG_SPECTEXT);
    printf("  Encoded: %u bytes (output: %u samples at 48 kHz)\n",
           enc_spec, spec_out_len);

    memset(recovered, 0, sizeof(recovered));
    unsigned dec_spec = MD_steg_decode(stego_spec, spec_out_len, 48000.0,
                                        recovered, sizeof(recovered),
                                        MD_STEG_SPECTEXT);
    printf("  Decoded: %u bytes -> \"%s\"\n", dec_spec, recovered);

    if (dec_spec != enc_spec || strcmp(recovered, spec_secret) != 0) {
        printf("  Spectext FAILED: decoded message does not match!\n");
        pass = 0;
    } else {
        printf("  Spectext PASSED\n");
    }
    free(stego_spec);

    if (pass)
        printf("\nSelf-test PASSED: all methods recovered the message.\n");
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
        "  %s                                              (self-test)\n"
        "  %s FILE                                         (auto-detect)\n"
        "  %s --encode METHOD MSG [-i HOST] [-o OUT]\n"
        "  %s --decode [METHOD] FILE\n"
        "  %s --encode-image METHOD IMAGE [-i HOST] [-o OUT]\n"
        "  %s --decode-image [METHOD] FILE -o IMAGE_OUT\n"
        "\n"
        "METHOD: \"lsb\", \"freq\", or \"spectext\" (optional for decode — auto-detects if omitted)\n"
        "\n"
        "Examples:\n"
        "  %s --encode lsb \"secret message\" -o stego.wav\n"
        "  %s --decode stego.wav\n"
        "  %s stego.wav\n"
        "  %s --encode-image lsb space_invader.png -o stego.wav\n"
        "  %s --decode-image stego.wav -o recovered.png\n",
        prog, prog, prog, prog, prog, prog, prog, prog, prog, prog, prog);
}

int main(int argc, char *argv[])
{
    if (argc == 1)
        return self_test();

    if (argc >= 4 && strcmp(argv[1], "--encode") == 0) {
        int method = parse_method(argv[2]);
        if (method < 0) {
            fprintf(stderr, "Unknown method '%s' (use 'lsb', 'freq', or 'spectext')\n",
                    argv[2]);
            return 1;
        }
        const char *message = argv[3];
        const char *infile  = NULL;
        const char *outfile = "steg_output.wav";

        for (int i = 4; i < argc - 1; i++) {
            if (strcmp(argv[i], "-i") == 0) infile  = argv[++i];
            else if (strcmp(argv[i], "-o") == 0) outfile = argv[++i];
        }
        return encode_wav(method, message, infile, outfile);
    }

    if (argc >= 3 && strcmp(argv[1], "--decode") == 0) {
        int method = parse_method(argv[2]);
        if (method >= 0) {
            /* --decode METHOD FILE */
            if (argc < 4) { usage(argv[0]); return 1; }
            return decode_wav(method, argv[3]);
        }
        /* --decode FILE (auto-detect) */
        return auto_decode_wav(argv[2]);
    }

    if (argc >= 4 && strcmp(argv[1], "--encode-image") == 0) {
        int method = parse_method(argv[2]);
        if (method < 0) {
            fprintf(stderr, "Unknown method '%s' (use 'lsb', 'freq', or 'spectext')\n",
                    argv[2]);
            return 1;
        }
        const char *image_path = argv[3];
        const char *infile  = NULL;
        const char *outfile = "steg_output.wav";

        for (int i = 4; i < argc - 1; i++) {
            if (strcmp(argv[i], "-i") == 0) infile  = argv[++i];
            else if (strcmp(argv[i], "-o") == 0) outfile = argv[++i];
        }
        return encode_image_wav(method, image_path, infile, outfile);
    }

    if (argc >= 3 && strcmp(argv[1], "--decode-image") == 0) {
        int method = parse_method(argv[2]);
        const char *infile;
        int opt_start;

        if (method >= 0) {
            /* --decode-image METHOD FILE -o OUT */
            if (argc < 4) { usage(argv[0]); return 1; }
            infile = argv[3];
            opt_start = 4;
        } else {
            /* --decode-image FILE -o OUT (auto-detect) */
            infile = argv[2];
            opt_start = 3;

            /* Read WAV to detect method. */
            double  *stego    = NULL;
            unsigned signal_len, samprate;
            if (read_wav_to_double(infile, &stego, &signal_len, &samprate) != 0)
                return 1;
            method = MD_steg_detect(stego, signal_len, (double)samprate,
                                    NULL);
            free(stego);
            if (method < 0) {
                fprintf(stderr, "No hidden payload detected in %s\n", infile);
                return 1;
            }
            printf("Auto-detected method: %s\n", method_name(method));
        }

        const char *outfile = NULL;
        for (int i = opt_start; i < argc - 1; i++) {
            if (strcmp(argv[i], "-o") == 0) outfile = argv[++i];
        }
        if (!outfile) {
            fprintf(stderr, "Error: --decode-image requires -o OUTPUT_FILE\n");
            return 1;
        }
        return decode_image_wav(method, infile, outfile);
    }

    /* Bare filename: ./audio_steg FILE */
    if (argc == 2 && argv[1][0] != '-')
        return auto_decode_wav(argv[1]);

    usage(argv[0]);
    return 1;
}
