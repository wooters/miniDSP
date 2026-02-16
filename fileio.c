/**
 * @file fileio.c
 * @brief Audio file I/O and HTK feature file writing.
 *
 * This module reads audio using libsndfile (which supports dozens of
 * formats including WAV, FLAC, AIFF, and Ogg) and writes feature
 * vectors in HTK binary format for speech-processing pipelines.
 */
#include "fileio.h"

/* -----------------------------------------------------------------------
 * HTK file format support
 *
 * HTK (Hidden Markov Model Toolkit) is a classic speech recognition
 * framework.  Its binary feature file format starts with a 12-byte
 * header, followed by the feature data.  All multi-byte values are
 * stored in big-endian byte order.
 * -----------------------------------------------------------------------*/

/**
 * HTK file header structure (12 bytes total).
 *
 * Fields:
 *   nvecs      - Number of feature vectors in the file.
 *   sampperiod - Sampling period in units of 100 nanoseconds.
 *   vecsize    - Size of each vector in *bytes* (not floats).
 *   parmkind   - Parameter kind code (9 = USER-defined).
 */
typedef struct FIO_HTKheader {
    uint32_t nvecs;
    uint32_t sampperiod;
    uint16_t vecsize;
    uint16_t parmkind;
} FIO_HTKheader;

/* -----------------------------------------------------------------------
 * Byte-order helpers
 *
 * HTK files are always big-endian.  Most modern PCs are little-endian,
 * so we need to swap bytes before writing.  These helpers detect the
 * host byte order at runtime and swap if necessary.
 * -----------------------------------------------------------------------*/

/** Check if the host is big-endian. */
static int is_bigendian(void)
{
    const int i = 1;
    return (*(const char *)&i) == 0;
}

/** Reverse the bytes of a value in-place (e.g. convert little-endian to big). */
static void swap_bytes(void *pv, size_t n)
{
    char *p = (char *)pv;
    for (size_t lo = 0, hi = n - 1; hi > lo; lo++, hi--) {
        char tmp = p[lo];
        p[lo] = p[hi];
        p[hi] = tmp;
    }
}

/** Convenience macro: swap the bytes of variable x. */
#define SWAP(x) swap_bytes(&(x), sizeof(x))

/* -----------------------------------------------------------------------
 * Internal helpers
 * -----------------------------------------------------------------------*/

/**
 * Allocate memory or terminate the program on failure.
 * This is a safety net -- in a small CLI tool it is better to crash
 * immediately with a clear message than to limp along with NULL pointers.
 */
static void *malloc_or_die(size_t nbytes, const char *msg)
{
    void *tmp = malloc(nbytes);
    if (tmp == NULL) {
        fprintf(stderr, "%s", msg);
        exit(1);
    }
    return tmp;
}

/* -----------------------------------------------------------------------
 * Public API
 * -----------------------------------------------------------------------*/

/**
 * Read a single-channel audio file into a float array.
 *
 * Supported formats include everything that libsndfile can open:
 * WAV, FLAC, AIFF, OGG, and many more.
 *
 * @param infile    Path to the audio file.
 * @param indata    Output: pointer to allocated float array with samples.
 *                  The caller is responsible for calling free() on this.
 * @param datalen   Output: total number of samples read.
 * @param samprate  Output: sample rate in Hz.
 * @param donorm    1 to normalise samples to [-1.0, 1.0]; 0 for raw values.
 */
void FIO_read_audio(const char *infile,
                    float **indata,
                    size_t *datalen,
                    unsigned *samprate,
                    unsigned donorm)
{
    SF_INFO sfinfo;
    memset(&sfinfo, 0, sizeof(sfinfo));

    /* Open the file and read its metadata */
    SNDFILE *sf = sf_open(infile, SFM_READ, &sfinfo);
    if (sf == NULL) {
        fprintf(stderr, "Error opening audio file: %s\n", infile);
        exit(1);
    }

    if (sfinfo.channels != 1) {
        fprintf(stderr,
                "Input file has %d channels. Only mono files are supported.\n"
                "Use 'sox' to split multi-channel files.\n",
                sfinfo.channels);
        sf_close(sf);
        exit(1);
    }

    sf_count_t nsamps = sfinfo.frames;

    /* Tell libsndfile whether to normalise float output to [-1, 1] */
    if (donorm == 1)
        sf_command(sf, SFC_SET_NORM_FLOAT, NULL, SF_TRUE);
    else
        sf_command(sf, SFC_SET_NORM_FLOAT, NULL, SF_FALSE);

    /* Allocate and read */
    float *tmpdata = malloc_or_die((size_t)nsamps * sizeof(float),
                                   "Error allocating memory for audio data\n");

    sf_count_t nread = sf_read_float(sf, tmpdata, nsamps);
    if (nread != nsamps) {
        fprintf(stderr, "Error reading %s: expected %ld samples, got %ld\n",
                infile, (long)nsamps, (long)nread);
        free(tmpdata);
        sf_close(sf);
        exit(1);
    }

    sf_close(sf);

    /* Set output parameters */
    *indata  = tmpdata;
    *datalen = (size_t)nsamps;
    *samprate = (unsigned)sfinfo.samplerate;
}

/**
 * Write feature vectors in HTK binary file format.
 *
 * The HTK format is:
 *   [12-byte header][vector 1][vector 2]...[vector N]
 *
 * Each vector is a sequence of big-endian 32-bit floats.
 *
 * @param outfile       Output file path.
 * @param outvecs       Array of nvecs pointers, each to veclen floats.
 * @param nvecs         Number of feature vectors.
 * @param veclen        Number of float elements per vector.
 * @param vecsamprate   Feature sampling rate in Hz.
 */
void FIO_write_htk_feats(const char *outfile,
                         const float **outvecs,
                         size_t nvecs,
                         size_t veclen,
                         unsigned vecsamprate)
{
    /* Build the 12-byte header */
    FIO_HTKheader hdr;
    hdr.nvecs      = (uint32_t)nvecs;
    hdr.sampperiod = (uint32_t)(1.0 / (float)vecsamprate * 1e7);
    hdr.vecsize    = (uint16_t)(veclen * sizeof(float));
    hdr.parmkind   = 9;  /* 9 = USER (user-defined parameter type) */

    FILE *f = fopen(outfile, "wb");
    if (f == NULL) {
        fprintf(stderr, "Error opening output file: %s\n", outfile);
        exit(1);
    }

    /* HTK files are big-endian, so swap bytes on little-endian machines */
    if (!is_bigendian()) {
        SWAP(hdr.nvecs);
        SWAP(hdr.sampperiod);
        SWAP(hdr.vecsize);
        SWAP(hdr.parmkind);
    }

    /* Write the header */
    fwrite(&hdr, sizeof(FIO_HTKheader), 1, f);

    /* Write the feature vectors */
    for (size_t i = 0; i < nvecs; i++) {
        if (!is_bigendian()) {
            /* Swap each float individually */
            for (size_t j = 0; j < veclen; j++) {
                float tmp = outvecs[i][j];
                SWAP(tmp);
                fwrite(&tmp, sizeof(float), 1, f);
            }
        } else {
            fwrite(outvecs[i], sizeof(float), veclen, f);
        }
    }

    fclose(f);
}
