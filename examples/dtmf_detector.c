/**
 * @file dtmf_detector.c
 * @brief Example: DTMF tone generation and detection.
 *
 * Three modes of operation:
 *
 *   1) No arguments — self-test: generate a DTMF sequence, detect it,
 *      and verify the results.
 *
 *   2) --generate DIGITS [-o FILE] — generate a DTMF tone sequence and
 *      write it to a WAV file.
 *
 *   3) --detect FILE — read a WAV file and detect DTMF tones.
 *
 * Build and run (from the repository root):
 *   make -C examples dtmf_detector
 *   cd examples && ./dtmf_detector
 *   cd examples && ./dtmf_detector --generate "5551234" -o tones.wav
 *   cd examples && ./dtmf_detector --detect tones.wav
 */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "minidsp.h"
#include "fileio.h"

/* --------------------------------------------------------------------- */

//! [generate-wav]
static int valid_dtmf_char(char ch)
{
    return (ch >= '0' && ch <= '9') || ch == '*' || ch == '#'
        || ch == 'A' || ch == 'a' || ch == 'B' || ch == 'b'
        || ch == 'C' || ch == 'c' || ch == 'D' || ch == 'd';
}

static int generate_wav(const char *digits, const char *outfile)
{
    const double sample_rate = 8000.0;
    const unsigned tone_ms   = 70;
    const unsigned pause_ms  = 70;

    for (const char *p = digits; *p; p++) {
        if (!valid_dtmf_char(*p)) {
            fprintf(stderr, "Invalid DTMF character '%c'. "
                            "Valid: 0-9, A-D, *, #\n", *p);
            return 1;
        }
    }

    unsigned num_digits = (unsigned)strlen(digits);
    unsigned signal_len = MD_dtmf_signal_length(num_digits, sample_rate,
                                                tone_ms, pause_ms);

    double *signal = malloc(signal_len * sizeof(double));
    if (!signal) { fprintf(stderr, "allocation failed\n"); return 1; }

    MD_dtmf_generate(signal, digits, sample_rate, tone_ms, pause_ms);

    /* Convert double -> float for WAV writing. */
    float *fdata = malloc(signal_len * sizeof(float));
    if (!fdata) { free(signal); fprintf(stderr, "allocation failed\n"); return 1; }
    for (unsigned i = 0; i < signal_len; i++)
        fdata[i] = (float)signal[i];

    int ret = FIO_write_wav(outfile, fdata, signal_len, (unsigned)sample_rate);
    if (ret == 0)
        printf("Generated DTMF \"%s\" -> %s  (%u samples, %.3f s)\n",
               digits, outfile, signal_len,
               (double)signal_len / sample_rate);
    else
        fprintf(stderr, "Error writing %s\n", outfile);

    free(fdata);
    free(signal);
    return ret;
}
//! [generate-wav]

/* --------------------------------------------------------------------- */

//! [detect-file]
static int detect_file(const char *infile)
{
    float  *fdata    = nullptr;
    size_t  datalen  = 0;
    unsigned samprate = 0;

    if (FIO_read_audio(infile, &fdata, &datalen, &samprate, 1) != 0) {
        fprintf(stderr, "Error reading %s\n", infile);
        return 1;
    }

    printf("Read %s: %zu samples at %u Hz (%.3f s)\n",
           infile, datalen, samprate, (double)datalen / (double)samprate);

    if (samprate < 4000) {
        fprintf(stderr, "Sample rate %u Hz is too low for DTMF detection "
                        "(minimum 4000 Hz)\n", samprate);
        free(fdata);
        return 1;
    }

    if (datalen > UINT_MAX) {
        fprintf(stderr, "File too large (%zu samples, max %u)\n",
                datalen, UINT_MAX);
        free(fdata);
        return 1;
    }

    /* Convert float -> double for the library. */
    double *signal = malloc(datalen * sizeof(double));
    if (!signal) {
        free(fdata);
        fprintf(stderr, "allocation failed\n");
        return 1;
    }
    for (size_t i = 0; i < datalen; i++)
        signal[i] = (double)fdata[i];
    free(fdata);

    /* Detect. */
    MD_DTMFTone tones[256];
    unsigned n = MD_dtmf_detect(signal, (unsigned)datalen,
                                (double)samprate, tones, 256);

    printf("\nDetected %u DTMF tone%s:\n", n, n == 1 ? "" : "s");
    if (n > 0) {
        printf("  %-6s  %-12s  %-12s\n", "Digit", "Start (s)", "End (s)");
        for (unsigned i = 0; i < n; i++)
            printf("  %-6c  %-12.3f  %-12.3f\n",
                   tones[i].digit, tones[i].start_s, tones[i].end_s);
    }

    free(signal);
    MD_shutdown();
    return 0;
}
//! [detect-file]

/* --------------------------------------------------------------------- */

//! [self-test]
static int self_test(void)
{
    const char    *test_digits = "14*258039#";
    const double   sample_rate = 8000.0;
    const unsigned tone_ms     = 70;
    const unsigned pause_ms    = 70;

    unsigned num_digits = (unsigned)strlen(test_digits);
    unsigned signal_len = MD_dtmf_signal_length(num_digits, sample_rate,
                                                tone_ms, pause_ms);

    printf("Self-test: generating DTMF sequence \"%s\"\n", test_digits);
    printf("  sample_rate = %.0f Hz, tone = %u ms, pause = %u ms\n",
           sample_rate, tone_ms, pause_ms);

    double *signal = malloc(signal_len * sizeof(double));
    if (!signal) { fprintf(stderr, "allocation failed\n"); return 1; }

    MD_dtmf_generate(signal, test_digits, sample_rate, tone_ms, pause_ms);

    MD_DTMFTone tones[64];
    unsigned n = MD_dtmf_detect(signal, signal_len, sample_rate, tones, 64);

    printf("\nDetected %u DTMF tone%s:\n", n, n == 1 ? "" : "s");
    printf("  %-6s  %-12s  %-12s\n", "Digit", "Start (s)", "End (s)");
    for (unsigned i = 0; i < n; i++)
        printf("  %-6c  %-12.3f  %-12.3f\n",
               tones[i].digit, tones[i].start_s, tones[i].end_s);

    /* Verify. */
    int pass = 1;
    if (n != num_digits) {
        printf("\nSelf-test FAILED: expected %u digits, detected %u\n",
               num_digits, n);
        pass = 0;
    } else {
        for (unsigned i = 0; i < num_digits; i++) {
            if (tones[i].digit != test_digits[i]) {
                printf("\nSelf-test FAILED: digit %u expected '%c' got '%c'\n",
                       i, test_digits[i], tones[i].digit);
                pass = 0;
                break;
            }
        }
    }

    if (pass)
        printf("\nSelf-test PASSED: all %u digits detected correctly\n",
               num_digits);

    free(signal);
    MD_shutdown();
    return pass ? 0 : 1;
}
//! [self-test]

/* --------------------------------------------------------------------- */

int main(int argc, char *argv[])
{
    if (argc == 1)
        return self_test();

    if (argc >= 3 && strcmp(argv[1], "--generate") == 0) {
        const char *digits  = argv[2];
        const char *outfile = "dtmf_output.wav";
        if (argc >= 5 && strcmp(argv[3], "-o") == 0)
            outfile = argv[4];
        return generate_wav(digits, outfile);
    }

    if (argc >= 3 && strcmp(argv[1], "--detect") == 0)
        return detect_file(argv[2]);

    fprintf(stderr,
            "Usage: %s\n"
            "       %s --generate DIGITS [-o FILE]\n"
            "       %s --detect FILE\n",
            argv[0], argv[0], argv[0]);
    return 1;
}
