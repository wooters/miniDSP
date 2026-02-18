/**
 * @file liveio.h
 * @brief Live audio recording and playback via PortAudio.
 *
 * This module wraps the PortAudio library to provide simple functions
 * for recording from a microphone and playing back to speakers.
 * All functions are non-blocking -- they return immediately while audio
 * I/O happens in the background via callbacks.
 *
 * Typical usage:
 *   1. Call LA_init() to initialise the audio system.
 *   2. Call LA_record() to start recording into a buffer.
 *   3. Poll LA_is_recording() until recording finishes.
 *   4. Call LA_play() to play back the recorded audio.
 *   5. Call LA_terminate() when done.
 */
#ifndef LIVEIO_H
#define LIVEIO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <portaudio.h>

/** Return type for live audio functions. */
typedef int LaError_t;

/** Error codes returned by live audio functions. */
enum LA_ERRORCODE {
    LA_OK    = 0,  /**< Operation succeeded.     */
    LA_ERROR = 1   /**< Something went wrong.    */
};

/** Recording modes for LA_record(). */
enum LA_RECORDTYPE {
    LA_REC_ONCE = 0, /**< Record until the buffer is full, then stop. */
    LA_REC_CONT      /**< Record continuously in a circular buffer.   */
};

/** Initialise the PortAudio system.  Must be called before any other LA_ function. */
LaError_t LA_init(void);

/** Shut down the PortAudio system.  Call when you are finished with audio. */
LaError_t LA_terminate(void);

/** Print the last PortAudio error message to the given stream. */
void LA_print_last_error(FILE *stream);

/** Check if the recording stream is active. Returns 1 if recording, 0 if not, -1 on error. */
int LA_is_recording(void);

/** Get the current time (in seconds) of the recording stream.  Returns -1.0 if inactive. */
double LA_get_cur_record_time(void);

/** Stop the current recording. */
void LA_stop_recording(void);

/**
 * Record audio from the default input device into a buffer.
 *
 * @param buffer   Pre-allocated buffer of 16-bit signed integers.
 * @param bufsize  Number of samples the buffer can hold.
 * @param samprate Desired sampling rate in Hz (e.g. 16000, 44100).
 * @param rectype  LA_REC_ONCE (fill buffer then stop) or LA_REC_CONT (circular).
 * @return         LA_OK on success, LA_ERROR on failure.
 *
 * @note Non-blocking: returns immediately while recording continues in the background.
 */
LaError_t LA_record(void *buffer, unsigned long bufsize,
                    unsigned samprate, int rectype);

/**
 * Record using a user-supplied PortAudio callback function.
 *
 * @param samprate  Sampling rate in Hz.
 * @param cb        Your PortAudio callback function.
 * @param cb_data   Data pointer passed to your callback.
 * @return          LA_OK on success, LA_ERROR on failure.
 */
LaError_t LA_record_callback(unsigned samprate, PaStreamCallback *cb, void *cb_data);

/** Check if the playback stream is active. Returns 1 if playing, 0 if not, -1 on error. */
int LA_is_playing(void);

/** Get the current time (in seconds) of the playback stream.  Returns -1.0 if inactive. */
double LA_get_cur_play_time(void);

/** Stop the current playback. */
void LA_stop_playing(void);

/**
 * Play audio from a buffer to the default output device.
 *
 * @param buffer   Buffer of 16-bit signed integers.
 * @param bufsize  Number of samples in the buffer.
 * @param samprate Sampling rate in Hz.
 * @return         LA_OK on success, LA_ERROR on failure.
 *
 * @note Non-blocking: returns immediately while playback continues in the background.
 *       Mono input is duplicated to both stereo channels.
 */
LaError_t LA_play(const void *buffer, unsigned long bufsize, unsigned samprate);

/**
 * Play using a user-supplied PortAudio callback function.
 *
 * @param samprate  Sampling rate in Hz.
 * @param cb        Your PortAudio callback function.
 * @param cb_data   Data pointer passed to your callback.
 * @return          LA_OK on success, LA_ERROR on failure.
 */
LaError_t LA_play_callback(unsigned samprate, PaStreamCallback *cb, void *cb_data);

#endif /* LIVEIO_H */
