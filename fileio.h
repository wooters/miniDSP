/**
 * @file fileio.h
 * @brief Audio file I/O and HTK feature file writing.
 *
 * This module provides:
 *   - Reading audio files (WAV, FLAC, AIFF, etc.) via libsndfile.
 *   - Writing feature vectors in HTK binary format, which is used by
 *     speech recognition toolkits like HTK and Kaldi.
 */
#ifndef FILEIO_H
#define FILEIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <sndfile.h>

/**
 * Read a single-channel audio file into memory.
 *
 * @param infile    Path to the audio file.
 * @param indata    Output: pointer to the audio samples (caller must free).
 * @param datalen   Output: number of samples read.
 * @param samprate  Output: sampling rate of the file (Hz).
 * @param donorm    If 1, normalise samples to [-1.0, 1.0].
 *
 * @note Only single-channel (mono) files are supported.
 *       Use the 'sox' utility to split multi-channel files first.
 */
void FIO_read_audio(const char *infile,
                    float **indata,
                    size_t *datalen,
                    unsigned *samprate,
                    unsigned donorm);

/**
 * Write feature vectors in HTK binary format.
 *
 * @param outfile       Path to the output file.
 * @param outvecs       Array of nvecs pointers, each pointing to veclen floats.
 * @param nvecs         Number of feature vectors.
 * @param veclen        Number of floats per vector.
 * @param vecsamprate   Sampling rate of the feature vectors (Hz).
 */
void FIO_write_htk_feats(const char *outfile,
                         const float **outvecs,
                         size_t nvecs,
                         size_t veclen,
                         unsigned vecsamprate);

#endif /* FILEIO_H */
