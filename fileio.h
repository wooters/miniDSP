/**
 * @file fileio.h
 * @brief Audio file I/O and feature vector writing.
 *
 * This module provides:
 *   - Reading audio files (WAV, FLAC, AIFF, etc.) via libsndfile.
 *   - Writing audio files (WAV float) via libsndfile.
 *   - Writing feature vectors in NumPy .npy format for Python interop.
 *   - Writing feature vectors in safetensors format for ML pipelines.
 *   - Writing feature vectors in HTK binary format (deprecated).
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
 * Write a 2D float32 array in NumPy .npy v1.0 format.
 *
 * Produces a file readable by numpy.load(). Data is stored as
 * little-endian float32, row-major (C order).
 *
 * @param outfile  Path to the output .npy file.
 * @param outvecs  Array of nvecs pointers, each pointing to veclen floats.
 * @param nvecs    Number of feature vectors (rows).
 * @param veclen   Number of floats per vector (columns).
 */
void FIO_write_npy(const char *outfile,
                   const float **outvecs,
                   size_t nvecs,
                   size_t veclen);

/**
 * Write a 2D float32 array in safetensors format.
 *
 * Produces a file readable by the safetensors Python library.
 * The tensor is stored under the key "features" as little-endian float32.
 *
 * @param outfile  Path to the output .safetensors file.
 * @param outvecs  Array of nvecs pointers, each pointing to veclen floats.
 * @param nvecs    Number of feature vectors (rows).
 * @param veclen   Number of floats per vector (columns).
 */
void FIO_write_safetensors(const char *outfile,
                           const float **outvecs,
                           size_t nvecs,
                           size_t veclen);

/**
 * Write mono float audio to a WAV file.
 *
 * Uses IEEE float format (SF_FORMAT_FLOAT) for lossless DSP round-trips.
 *
 * @param outfile   Path to the output .wav file.
 * @param data      Audio samples.
 * @param datalen   Number of samples.
 * @param samprate  Sampling rate in Hz.
 */
void FIO_write_wav(const char *outfile,
                   const float *data,
                   size_t datalen,
                   unsigned samprate);

/**
 * Write feature vectors in HTK binary format.
 *
 * @deprecated Use FIO_write_npy() or FIO_write_safetensors() instead.
 *
 * @param outfile       Path to the output file.
 * @param outvecs       Array of nvecs pointers, each pointing to veclen floats.
 * @param nvecs         Number of feature vectors.
 * @param veclen        Number of floats per vector.
 * @param vecsamprate   Sampling rate of the feature vectors (Hz).
 */
[[deprecated("use FIO_write_npy or FIO_write_safetensors")]]
void FIO_write_htk_feats(const char *outfile,
                         const float **outvecs,
                         size_t nvecs,
                         size_t veclen,
                         unsigned vecsamprate);

#endif /* FILEIO_H */
