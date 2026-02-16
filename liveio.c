/**
 * @file liveio.c
 * @brief Live audio recording and playback implementation.
 * @author Chuck Wooters <wooters@hey.com>
 *
 * This is a wrapper around PortAudio (http://portaudio.com), a
 * cross-platform library for real-time audio I/O.
 *
 * PortAudio works on a callback model: you give it a function, and it
 * calls that function whenever it needs more audio data (for playback)
 * or has new data available (for recording).  The callbacks run on a
 * separate high-priority thread managed by the OS audio subsystem.
 */

#include "liveio.h"

/* -----------------------------------------------------------------------
 * Internal data structures
 * -----------------------------------------------------------------------*/

/**
 * Simple audio buffer with read/write positions.
 *
 * This struct is shared between the main thread and the PortAudio
 * callback thread.  The callback writes (during recording) or reads
 * (during playback) from the buffer.
 */
typedef struct {
    int16_t      *audiodata;  /**< Pointer to the sample buffer.          */
    unsigned long len;        /**< Total number of samples in the buffer. */
    unsigned long pos_r;      /**< Current read position (playback).      */
    unsigned long pos_w;      /**< Current write position (recording).    */
} LA_simpleData;

/* -----------------------------------------------------------------------
 * Module-level state
 * -----------------------------------------------------------------------*/

static int           LA_InitDone     = 0;         /**< Has PA been initialised?  */
static PaError       LA_LastError    = paNoError;  /**< Most recent PA error code */
static PaStream     *LA_PlayStream   = NULL;       /**< Active playback stream    */
static PaStream     *LA_RecordStream = NULL;       /**< Active recording stream   */
static LA_simpleData LA_data;                      /**< Shared buffer descriptor  */

/* -----------------------------------------------------------------------
 * Stream status queries
 * -----------------------------------------------------------------------*/

double LA_get_cur_record_time(void)
{
    if (Pa_IsStreamActive(LA_RecordStream) == 1) {
        return Pa_GetStreamTime(LA_RecordStream);
    }
    return -1.0;
}

double LA_get_cur_play_time(void)
{
    if (Pa_IsStreamActive(LA_PlayStream) == 1) {
        return Pa_GetStreamTime(LA_PlayStream);
    }
    return -1.0;
}

void LA_stop_recording(void)
{
    if (LA_is_recording() == 1) {
        Pa_CloseStream(LA_RecordStream);
    }
}

void LA_stop_playing(void)
{
    if (LA_is_playing() == 1) {
        Pa_CloseStream(LA_PlayStream);
    }
}

/**
 * Check if the recording stream is currently active.
 * @return 1 if recording, 0 if not, -1 on error.
 */
int LA_is_recording(void)
{
    int res = Pa_IsStreamActive(LA_RecordStream);
    if (res == 1 || res == 0) return res;
    return -1;
}

/**
 * Check if the playback stream is currently active.
 * @return 1 if playing, 0 if not, -1 on error.
 */
int LA_is_playing(void)
{
    int res = Pa_IsStreamActive(LA_PlayStream);
    if (res == 1 || res == 0) return res;
    return -1;
}

/* -----------------------------------------------------------------------
 * PortAudio callbacks
 *
 * These functions are called by PortAudio on a dedicated audio thread.
 * They must be fast and must not block (no malloc, no printf, no locks).
 * -----------------------------------------------------------------------*/

/**
 * Callback for continuous (circular buffer) recording.
 *
 * When the write position reaches the end of the buffer, it wraps
 * around to the beginning.  This means old data is overwritten --
 * the buffer always contains the most recent audio.
 */
static int _ContinuousRecordCallback(const void *inputBuffer,
                                     void *outputBuffer,
                                     unsigned long framesPerBuffer,
                                     const PaStreamCallbackTimeInfo *timeInfo,
                                     PaStreamCallbackFlags statusFlags,
                                     void *userData)
{
    LA_simpleData *data = (LA_simpleData *)userData;
    const int16_t *rptr = (const int16_t *)inputBuffer;

    /* Suppress unused-parameter warnings */
    (void)outputBuffer;
    (void)timeInfo;
    (void)statusFlags;

    for (unsigned long i = 0; i < framesPerBuffer; i++) {
        data->audiodata[data->pos_w] = (inputBuffer != NULL) ? *rptr++ : 0;
        data->pos_w++;
        if (data->pos_w >= data->len) {
            data->pos_w = 0;  /* Wrap around */
        }
    }

    return paContinue;
}

/**
 * Callback for one-shot recording (fill the buffer once, then stop).
 *
 * Returns paComplete when the buffer is full, which tells PortAudio
 * to stop the stream after this callback returns.
 */
static int _RecordCallback(const void *inputBuffer,
                           void *outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo *timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void *userData)
{
    LA_simpleData *data = (LA_simpleData *)userData;
    const int16_t *rptr = (const int16_t *)inputBuffer;
    int16_t       *wptr = &data->audiodata[data->pos_w];

    (void)outputBuffer;
    (void)timeInfo;
    (void)statusFlags;

    /* How many frames can we still fit in the buffer? */
    unsigned long framesLeft  = data->len - data->pos_w;
    unsigned long framesToRec;
    int finished;

    if (framesLeft <= framesPerBuffer) {
        framesToRec = framesLeft;
        finished    = paComplete;
    } else {
        framesToRec = framesPerBuffer;
        finished    = paContinue;
    }

    /* Copy samples (or write silence if input is NULL) */
    for (unsigned long i = 0; i < framesToRec; i++) {
        *wptr++ = (inputBuffer != NULL) ? *rptr++ : 0;
    }
    data->pos_w += framesToRec;

    return finished;
}

/**
 * Callback for playback: reads mono samples and duplicates them
 * to both stereo output channels (left and right get the same data).
 */
static int _PlayCallback(const void *inputBuffer,
                         void *outputBuffer,
                         unsigned long framesPerBuffer,
                         const PaStreamCallbackTimeInfo *timeInfo,
                         PaStreamCallbackFlags statusFlags,
                         void *userData)
{
    LA_simpleData *data = (LA_simpleData *)userData;
    int16_t       *out  = (int16_t *)outputBuffer;

    (void)timeInfo;
    (void)statusFlags;
    (void)inputBuffer;

    int finished = paContinue;
    unsigned long nout = framesPerBuffer;

    /* Check if we are near the end of the buffer */
    if (data->pos_r + framesPerBuffer >= data->len) {
        finished = paComplete;
        nout = data->len - data->pos_r;
    }

    /* Write each mono sample to both stereo channels */
    for (unsigned long i = 0; i < nout; i++) {
        *out++ = data->audiodata[data->pos_r];  /* Left channel  */
        *out++ = data->audiodata[data->pos_r];  /* Right channel */
        data->pos_r++;
    }

    return finished;
}

/**
 * Called by PortAudio when a one-shot recording stream finishes.
 * Closes the stream so resources are released.
 */
static void _RecordStreamFinished(void *userData)
{
    (void)userData;
    PaError err = Pa_CloseStream(LA_RecordStream);
    if (err != paNoError) {
        fprintf(stderr, "PortAudio Error: %s\n", Pa_GetErrorText(err));
    }
}

/**
 * Called by PortAudio when a playback stream finishes.
 * Closes the stream so resources are released.
 */
static void _PlayStreamFinished(void *userData)
{
    (void)userData;
    PaError err = Pa_CloseStream(LA_PlayStream);
    if (err != paNoError) {
        fprintf(stderr, "PortAudio Error: %s\n", Pa_GetErrorText(err));
    }
}

/* -----------------------------------------------------------------------
 * Public API: recording
 * -----------------------------------------------------------------------*/

LaError_t LA_record_callback(unsigned samprate, PaStreamCallback *cb, void *cb_data)
{
    /* Auto-initialise if needed */
    if (LA_InitDone == 0) {
        if (LA_init() != LA_OK) return LA_ERROR;
    }
    if (LA_is_recording() == 1) return LA_ERROR;

    PaStreamParameters inputParams;
    inputParams.device = Pa_GetDefaultInputDevice();
    if (inputParams.device == paNoDevice) {
        LA_LastError = paNoDevice;
        return LA_ERROR;
    }
    inputParams.channelCount            = 1;       /* Mono input */
    inputParams.sampleFormat            = paInt16;
    inputParams.suggestedLatency        = Pa_GetDeviceInfo(inputParams.device)->defaultHighInputLatency;
    inputParams.hostApiSpecificStreamInfo = NULL;

    PaError err;
    err = Pa_OpenStream(&LA_RecordStream, &inputParams, NULL,
                        samprate, paFramesPerBufferUnspecified,
                        paClipOff, cb, cb_data);
    if (err != paNoError) return LA_ERROR;

    err = Pa_StartStream(LA_RecordStream);
    if (err != paNoError) return LA_ERROR;

    return LA_OK;
}

LaError_t LA_record(void *buffer, unsigned long bufsize,
                    unsigned samprate, int rectype)
{
    if (LA_InitDone == 0) {
        if (LA_init() != LA_OK) return LA_ERROR;
    }
    if (LA_is_recording() == 1) return LA_ERROR;

    LA_data.audiodata = (int16_t *)buffer;
    LA_data.len       = bufsize;
    LA_data.pos_w     = 0;

    PaStreamParameters inputParams;
    inputParams.device = Pa_GetDefaultInputDevice();
    if (inputParams.device == paNoDevice) {
        LA_LastError = paNoDevice;
        return LA_ERROR;
    }
    inputParams.channelCount            = 1;
    inputParams.sampleFormat            = paInt16;
    inputParams.suggestedLatency        = Pa_GetDeviceInfo(inputParams.device)->defaultHighInputLatency;
    inputParams.hostApiSpecificStreamInfo = NULL;

    PaError err;
    PaStreamCallback *cb;

    if (rectype == LA_REC_ONCE) {
        cb = _RecordCallback;
    } else if (rectype == LA_REC_CONT) {
        cb = _ContinuousRecordCallback;
    } else {
        return LA_ERROR;
    }

    err = Pa_OpenStream(&LA_RecordStream, &inputParams, NULL,
                        samprate, paFramesPerBufferUnspecified,
                        paClipOff, cb, &LA_data);
    if (err != paNoError) return LA_ERROR;

    /* For one-shot recording, register a finished callback to auto-close */
    if (rectype == LA_REC_ONCE) {
        err = Pa_SetStreamFinishedCallback(LA_RecordStream, _RecordStreamFinished);
        if (err != paNoError) return LA_ERROR;
    }

    err = Pa_StartStream(LA_RecordStream);
    if (err != paNoError) return LA_ERROR;

    return LA_OK;
}

/* -----------------------------------------------------------------------
 * Public API: playback
 * -----------------------------------------------------------------------*/

LaError_t LA_play_callback(unsigned samprate, PaStreamCallback *cb, void *cb_data)
{
    if (LA_InitDone == 0) {
        if (LA_init() != LA_OK) return LA_ERROR;
    }
    if (LA_is_playing() == 1) return LA_ERROR;

    PaStreamParameters outputParams;
    outputParams.device = Pa_GetDefaultOutputDevice();
    if (outputParams.device == paNoDevice) {
        LA_LastError = paNoDevice;
        return LA_ERROR;
    }
    outputParams.channelCount            = 2;       /* Stereo output */
    outputParams.sampleFormat            = paInt16;
    outputParams.suggestedLatency        = Pa_GetDeviceInfo(outputParams.device)->defaultHighOutputLatency;
    outputParams.hostApiSpecificStreamInfo = NULL;

    PaError err;
    err = Pa_OpenStream(&LA_PlayStream, NULL, &outputParams,
                        samprate, paFramesPerBufferUnspecified,
                        paNoFlag, cb, cb_data);
    if (err != paNoError) return LA_ERROR;

    err = Pa_StartStream(LA_PlayStream);
    if (err != paNoError) return LA_ERROR;

    return LA_OK;
}

LaError_t LA_play(const void *buffer, unsigned long bufsize, unsigned samprate)
{
    if (LA_InitDone == 0) {
        if (LA_init() != LA_OK) return LA_ERROR;
    }
    if (LA_is_playing() == 1) return LA_ERROR;

    LA_data.audiodata = (int16_t *)buffer;
    LA_data.len       = bufsize;
    LA_data.pos_r     = 0;

    PaStreamParameters outputParams;
    outputParams.device = Pa_GetDefaultOutputDevice();
    if (outputParams.device == paNoDevice) {
        LA_LastError = paNoDevice;
        return LA_ERROR;
    }
    outputParams.channelCount            = 2;       /* Stereo output */
    outputParams.sampleFormat            = paInt16;
    outputParams.suggestedLatency        = Pa_GetDeviceInfo(outputParams.device)->defaultHighOutputLatency;
    outputParams.hostApiSpecificStreamInfo = NULL;

    PaError err;
    err = Pa_OpenStream(&LA_PlayStream, NULL, &outputParams,
                        samprate, paFramesPerBufferUnspecified,
                        paNoFlag, _PlayCallback, &LA_data);
    if (err != paNoError) return LA_ERROR;

    err = Pa_SetStreamFinishedCallback(LA_PlayStream, _PlayStreamFinished);
    if (err != paNoError) return LA_ERROR;

    err = Pa_StartStream(LA_PlayStream);
    if (err != paNoError) return LA_ERROR;

    return LA_OK;
}

/* -----------------------------------------------------------------------
 * Public API: error reporting and lifecycle
 * -----------------------------------------------------------------------*/

void LA_print_last_error(FILE *stream)
{
    fprintf(stream, "PortAudio Error: %s\n", Pa_GetErrorText(LA_LastError));
}

/**
 * Initialise the PortAudio system.
 *
 * If already initialised, terminates first and re-initialises.
 * @return LA_OK on success, LA_ERROR on failure.
 */
LaError_t LA_init(void)
{
    if (LA_InitDone) {
        if (LA_terminate() != LA_OK) return LA_ERROR;
    }

    LA_LastError = Pa_Initialize();
    if (LA_LastError != paNoError) {
        LA_InitDone = 0;
        return LA_ERROR;
    }

    memset(&LA_data, 0, sizeof(LA_data));
    LA_InitDone = 1;
    return LA_OK;
}

/**
 * Terminate the PortAudio system and release resources.
 * @return LA_OK on success, LA_ERROR on failure.
 */
LaError_t LA_terminate(void)
{
    if (LA_InitDone == 0) {
        LA_LastError = paNotInitialized;
        return LA_ERROR;
    }

    LA_LastError = Pa_Terminate();
    if (LA_LastError != paNoError) {
        LA_InitDone = 0;
        return LA_ERROR;
    }

    /* BUG FIX: previously this was incorrectly set to 1 after terminating,
     * which meant the system thought it was still initialised. */
    LA_InitDone = 0;
    return LA_OK;
}
