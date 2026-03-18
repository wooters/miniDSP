/**
 * @file resample.c
 * @brief Command-line sample rate conversion tool.
 *
 * Usage: resample [-z zero_crossings] [-b kaiser_beta] <input> <rate> <output.wav>
 *
 * Reads a mono audio file, resamples to the target rate using polyphase sinc
 * interpolation, and writes a WAV output file (IEEE float format).
 */

#include "minidsp.h"
#include "fileio.h"

#define DEFAULT_ZERO_CROSSINGS 32
#define DEFAULT_KAISER_BETA    10.0

static void usage(const char *prog)
{
    fprintf(stderr,
            "Usage: %s [-z zero_crossings] [-b kaiser_beta] <input> <rate> <output.wav>\n"
            "\n"
            "  Resample a mono audio file to a new sample rate.\n"
            "\n"
            "  <input>   Input audio file (mono; WAV, FLAC, AIFF, OGG, etc.)\n"
            "  <rate>    Target sample rate in Hz (positive integer)\n"
            "  <output>  Output file path (WAV format only)\n"
            "\n"
            "Options:\n"
            "  -z N      Number of sinc zero-crossings per side (default: %d)\n"
            "  -b F      Kaiser window beta parameter (default: %.1f)\n",
            prog, DEFAULT_ZERO_CROSSINGS, DEFAULT_KAISER_BETA);
}

/** Return 1 if str ends with ".wav" (case-insensitive). */
static int ends_with_wav(const char *str)
{
    size_t len = strlen(str);
    if (len < 4) return 0;
    const char *ext = str + len - 4;
    return (ext[0] == '.' &&
            (ext[1] == 'w' || ext[1] == 'W') &&
            (ext[2] == 'a' || ext[2] == 'A') &&
            (ext[3] == 'v' || ext[3] == 'V'));
}

int main(int argc, char *argv[])
{
    unsigned zero_crossings = DEFAULT_ZERO_CROSSINGS;
    double kaiser_beta = DEFAULT_KAISER_BETA;

    /* Parse optional flags */
    int i = 1;
    while (i < argc && argv[i][0] == '-') {
        if (strcmp(argv[i], "-z") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 1; }
            zero_crossings = (unsigned)atoi(argv[++i]);
            if (zero_crossings == 0) {
                fprintf(stderr, "Error: zero-crossings must be a positive integer\n");
                return 1;
            }
        } else if (strcmp(argv[i], "-b") == 0) {
            if (i + 1 >= argc) { usage(argv[0]); return 1; }
            kaiser_beta = atof(argv[++i]);
            if (kaiser_beta <= 0.0) {
                fprintf(stderr, "Error: kaiser beta must be positive\n");
                return 1;
            }
        } else {
            fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
        i++;
    }

    /* Remaining args: input, rate, output */
    if (argc - i != 3) {
        usage(argv[0]);
        return 1;
    }

    const char *input_path  = argv[i];
    int target_rate = atoi(argv[i + 1]);
    const char *output_path = argv[i + 2];

    if (target_rate <= 0) {
        fprintf(stderr, "Error: target sample rate must be a positive integer\n");
        return 1;
    }

    if (!ends_with_wav(output_path)) {
        fprintf(stderr, "Warning: output will be WAV format regardless of extension\n");
    }

    /* Read input audio */
    float *in_float = NULL;
    size_t in_len = 0;
    unsigned in_rate = 0;

    if (FIO_read_audio(input_path, &in_float, &in_len, &in_rate, 1) != 0) {
        return 1;
    }

    double in_duration = (double)in_len / (double)in_rate;

    /* Same-rate detection */
    if ((unsigned)target_rate == in_rate) {
        fprintf(stdout,
                "Note: input is already at %u Hz, copying without resampling\n",
                in_rate);

        if (FIO_write_wav(output_path, in_float, in_len, in_rate) != 0) {
            free(in_float);
            return 1;
        }

        printf("Input:  %s (%u Hz, %zu samples, %.2f s)\n",
               input_path, in_rate, in_len, in_duration);
        printf("Output: %s (%u Hz, %zu samples)\n",
               output_path, in_rate, in_len);

        free(in_float);
        return 0;
    }

    /* Convert float -> double */
    double *in_double = malloc(in_len * sizeof(double));
    if (in_double == NULL) {
        fprintf(stderr, "Error: failed to allocate memory\n");
        free(in_float);
        return 1;
    }
    for (size_t k = 0; k < in_len; k++) {
        in_double[k] = (double)in_float[k];
    }
    free(in_float);

    //! [resample-core]
    /* Resample */
    unsigned out_len = MD_resample_output_len(
        (unsigned)in_len, (double)in_rate, (double)target_rate);

    double *out_double = malloc(out_len * sizeof(double));
    if (out_double == NULL) {
        fprintf(stderr, "Error: failed to allocate memory\n");
        free(in_double);
        return 1;
    }

    MD_resample(in_double, (unsigned)in_len,
                out_double, out_len,
                (double)in_rate, (double)target_rate,
                zero_crossings, kaiser_beta);
    //! [resample-core]
    free(in_double);

    /* Convert double -> float */
    float *out_float = malloc(out_len * sizeof(float));
    if (out_float == NULL) {
        fprintf(stderr, "Error: failed to allocate memory\n");
        free(out_double);
        return 1;
    }
    for (unsigned k = 0; k < out_len; k++) {
        out_float[k] = (float)out_double[k];
    }
    free(out_double);

    /* Write output WAV */
    if (FIO_write_wav(output_path, out_float, out_len, (unsigned)target_rate) != 0) {
        free(out_float);
        return 1;
    }
    free(out_float);

    /* Print summary */
    printf("Input:  %s (%u Hz, %zu samples, %.2f s)\n",
           input_path, in_rate, in_len, in_duration);
    printf("Output: %s (%d Hz, %u samples)\n",
           output_path, target_rate, out_len);

    return 0;
}
