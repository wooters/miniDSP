/**
 * @file liveio.c
 * Implementation of live audio play/record system. This is basically a wrapper
 *  the PortAudio library (http://portaudio.com)
 *
 * @author Chuck Wooters <wooters@icsi.berkeley.edu>
 */

#include "liveio.h"

typedef struct
{
  int16_t* audiodata;   /**< buffer pointer */
  unsigned long len;    /**< length of buffer pointed to by audiodata */
  unsigned long pos_r;  /**< current read position (for playing) */
  unsigned long pos_w;  /**< current write position (for recording) */
} LA_simpleData;

static int     LA_InitDone  = 0;          /**< Has local audio been initialized? */
static PaError LA_LastError = paNoError;  /**< Error code for the last error encountered */
static PaStream* LA_PlayStream = NULL;    /**< Currently active play stream */
static PaStream* LA_RecordStream = NULL;  /**< Currently active record stream */
static LA_simpleData LA_data;             /**< Communication structure */

/**
 * Return the current time (in seconds) according the record stream.
 */
double LA_get_cur_record_time(void) {
  if (Pa_IsStreamActive(LA_RecordStream) == 1) {
    return Pa_GetStreamTime(LA_RecordStream);
  }
  return -1.0;
}

double LA_get_cur_play_time(void) {
  if (Pa_IsStreamActive(LA_PlayStream) == 1) {
    return Pa_GetStreamTime(LA_PlayStream);
  }
  return -1.0;
}

void LA_stop_recording(void) {
  if (LA_is_recording() == 1) {
    Pa_CloseStream(LA_RecordStream);
  }
}

void LA_stop_playing(void) {
  if (LA_is_playing() == 1) {
    Pa_CloseStream(LA_PlayStream);
  }
}

/**
 * A callback function for continuously recording audio into a buffer.
 */
static int _ContinuousRecordCallback(const void* inputBuffer, 
				     void* outputBuffer,
				     unsigned long framesPerBuffer,
				     const PaStreamCallbackTimeInfo* timeInfo,
				     PaStreamCallbackFlags statusFlags,
				     void* userData )
{
  LA_simpleData *data = (LA_simpleData*)userData;
  const int16_t* rptr = (const int16_t*)inputBuffer; /* read pointer */
  unsigned long i;

  (void) outputBuffer; /* Prevent unused variable warnings. */
  (void) timeInfo;
  (void) statusFlags;

  if( inputBuffer == NULL ) {
    for( i=0; i<framesPerBuffer; i++ ) {
      data->audiodata[data->pos_w++] = 0;
      if (data->pos_w >= data->len) data->pos_w=0;
    }
  } else {
    for( i=0; i<framesPerBuffer; i++ ) {
      data->audiodata[data->pos_w++] = *rptr++;
      if (data->pos_w >= data->len) data->pos_w=0;
    }
  }

  return paContinue;
}

/**
 * A callback function for recording audio into a buffer.
 */
static int _RecordCallback(const void* inputBuffer, 
			   void* outputBuffer,
                           unsigned long framesPerBuffer,
                           const PaStreamCallbackTimeInfo* timeInfo,
                           PaStreamCallbackFlags statusFlags,
                           void* userData )
{
  LA_simpleData *data = (LA_simpleData*)userData;
  const int16_t* rptr = (const int16_t*)inputBuffer;
  int16_t* wptr = &data->audiodata[data->pos_w];
  unsigned long framesToRec;
  unsigned long i;
  int finished;
  unsigned long framesLeft = data->len - data->pos_w - 1;

  (void) outputBuffer; /* Prevent unused variable warnings. */
  (void) timeInfo;
  (void) statusFlags;

  if( framesLeft < framesPerBuffer ) {
    framesToRec = framesLeft;
    finished = paComplete;
  } else {
    framesToRec = framesPerBuffer;
    finished = paContinue;
  }

  if( inputBuffer == NULL ) {
    for( i=0; i<framesToRec; i++ ) {
      *wptr++ = 0;
    }
  } else {
    for( i=0; i<framesToRec; i++ ) {
      *wptr++ = *rptr++;
    }
  }
  data->pos_w += framesToRec;
  return finished;
}


/**
 * A callback function for playing audio from a buffer.
 */
static int _PlayCallback( const void* inputBuffer, 
			  void* outputBuffer,
			  unsigned long framesPerBuffer,
			  const PaStreamCallbackTimeInfo* timeInfo,
			  PaStreamCallbackFlags statusFlags,
			  void* userData )
{
  LA_simpleData* data = (LA_simpleData*)userData;
  int16_t* out = (int16_t*) outputBuffer;
  (void) timeInfo;
  (void) statusFlags;
  (void) inputBuffer;

  int returnval=paContinue;
  unsigned long nout = framesPerBuffer;
  if ( (data->pos_r + framesPerBuffer) >= data->len) {
    returnval = paComplete;
    nout = data->len - data->pos_r -1;
  }

  for(unsigned long i=0; i<nout; i++) {
    *out++ = data->audiodata[data->pos_r];
    *out++ = data->audiodata[data->pos_r];
    data->pos_r += 1;
  }

  return returnval;
}

/**
 * A callback routine that is called when a record stream is
 * finished.
 */
static void _RecordStreamFinished( void* userData )
{
  int err = Pa_CloseStream(LA_RecordStream);
  if (err != paNoError) 
    printf("PortAudio Error: %s\n",Pa_GetErrorText(err));
}

/**
 * A callback routine that is called when a playback stream is
 * finished.
 */
static void _PlayStreamFinished( void* userData )
{
  int err = Pa_CloseStream(LA_PlayStream);
  if (err != paNoError) 
    printf("PortAudio Error: %s\n",Pa_GetErrorText(err));
}

/**
 * Determine if the recording stream is currently recording.
 *
 * @returns 1 if currently recording, 0 if not, and -1 if there was an error.
 */
int LA_is_recording(void) {
  int res = Pa_IsStreamActive(LA_RecordStream);
  if (res == 1 || res == 0) return res;
  return -1; /* indicates and error */
}

/**
 * Determine if the playback stream is currently playing something.
 *
 * @returns 1 if currently playing, 0 if not, and -1 if there was an error.
 */
int LA_is_playing(void) {
  int res = Pa_IsStreamActive(LA_PlayStream);
  if (res == 1 || res == 0) return res;
  return -1; /* indicates and error */
}

/**
 * Record from the mic using the user-specified callback and data structure.
 *
 * @param samprate Sampling rate of the audio
 * @param cb A callback function to pass to Pa_OpenStream()
 * @param cb_data A pointer to a data structure to be passed to the callback function
 *
 * @returns A value of type ::LaError_t where LA_OK indicates there were no problems, and a value of LA_ERROR
 * indicates that some error occurred. Use ::LA_print_last_error() to print the error.
 * 
 * @note This function is non-blocking and will return immediately (i.e. while it is still
 * recording. So use ::LA_is_recording() to determine when the sound is finished playing.
 *
 */
LaError_t LA_record_callback(const unsigned samprate, PaStreamCallback* cb, void* cb_data) {

  if (LA_InitDone == 0) { /* initialize portaudio system, if needed */
    int failed = LA_init();
    if (failed) return LA_ERROR;
  }
  if (LA_is_recording() == 1) return LA_ERROR;

  PaStreamParameters inputParameters;
  inputParameters.device = Pa_GetDefaultInputDevice();
  if (inputParameters.device == paNoDevice) {
    LA_LastError = paNoDevice;
    return LA_ERROR;
  }

  inputParameters.channelCount = 1;
  inputParameters.sampleFormat = paInt16; /* 16 bit ints */
  inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;

  PaError err;
  err = Pa_OpenStream(&LA_RecordStream,
		      &inputParameters,
		      NULL,
		      samprate,
		      paFramesPerBufferUnspecified,
		      paClipOff,
		      cb,
		      cb_data );

  if( err != paNoError ) return LA_ERROR;

  err = Pa_StartStream(LA_RecordStream);
  if (err != paNoError) return LA_ERROR;

  return LA_OK;
}

/**
 * Record from the mic and save into a buffer.
 *
 * @param buffer A buffer of 16-bit ints
 * @param bufsize Size of the \a buffer
 * @param samprate Sampling rate of the audio
 *
 * @returns A value of type ::LaError_t where LA_OK indicates there were no problems, and a value of LA_ERROR
 * indicates that some error occurred. Use ::LA_print_last_error() to print the error.
 * 
 * @note This function is non-blocking and will return immediately (i.e. while it is still
 * recording. So use ::LA_is_recording() to determine when the sound is finished playing.
 *
 */
LaError_t LA_record(void* const buffer, unsigned long bufsize, unsigned samprate, int rectype) {

  if (LA_InitDone == 0) { /* initialize portaudio system, if needed */
    int failed = LA_init();
    if (failed) return LA_ERROR;
  }
  if (LA_is_recording() == 1) return LA_ERROR;

  LA_data.audiodata = (int16_t*) buffer;
  LA_data.len = bufsize;
  LA_data.pos_w = 0;

  PaStreamParameters inputParameters;
  inputParameters.device = Pa_GetDefaultInputDevice();
  if (inputParameters.device == paNoDevice) {
    LA_LastError = paNoDevice;
    return LA_ERROR;
  }
  //printf( "Input device # %d.\n", inputParameters.device );
  //printf( "Input LL: %g s\n", Pa_GetDeviceInfo( inputParameters.device )->defaultLowInputLatency );
  //printf( "Input HL: %g s\n", Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency );
  inputParameters.channelCount = 1;
  inputParameters.sampleFormat = paInt16; /* 16 bit ints */
  inputParameters.suggestedLatency = Pa_GetDeviceInfo( inputParameters.device )->defaultHighInputLatency;
  inputParameters.hostApiSpecificStreamInfo = NULL;

  PaError err;
  if (rectype == LA_REC_ONCE) {
    err = Pa_OpenStream(&LA_RecordStream,
			&inputParameters,
			NULL,
			samprate,
			paFramesPerBufferUnspecified,
			paClipOff,
			_RecordCallback,
			&LA_data );
  } else if (rectype == LA_REC_CONT) {
    err = Pa_OpenStream(&LA_RecordStream,
			&inputParameters,
			NULL,
			samprate,
			paFramesPerBufferUnspecified,
			paClipOff,
			_ContinuousRecordCallback,
			&LA_data );
  } else {
    return LA_ERROR;
  }
  if( err != paNoError ) return LA_ERROR;

  if (rectype == LA_REC_ONCE) {
    err = Pa_SetStreamFinishedCallback(LA_RecordStream, _RecordStreamFinished);
    if (err != paNoError) return LA_ERROR;
  }

  err = Pa_StartStream(LA_RecordStream);
  if (err != paNoError) return LA_ERROR;

  return LA_OK;
}

/**
 * Play sound to the speakers using the user-specified callback and data structure.
 *
 * @param samprate Sampling rate of the audio
 * @param cb A callback function to pass to Pa_OpenStream()
 * @param cb_data A pointer to a data structure to be passed to the callback function
 *
 * @returns A value of type ::LaError_t where LA_OK indicates there were no problems, and a value of LA_ERROR
 * indicates that some error occurred. Use ::LA_print_last_error() to print the error.
 * 
 * @note This function is non-blocking and will return immediately (i.e. before the sound has 
 * finished playing. So use ::LA_is_playing() to determine when the sound is finished playing.
 *
 */
LaError_t LA_play_callback(const unsigned samprate, PaStreamCallback* cb, void* cb_data) {

  if (LA_InitDone == 0) { /* initialize portaudio system, if needed */
    int failed = LA_init();
    if (failed) return LA_ERROR;
  }
  if (LA_is_playing() == 1) return LA_ERROR;

  PaStreamParameters outputParameters;
  outputParameters.device = Pa_GetDefaultOutputDevice();
  //outputParameters.device = 4; // Soundflower (2ch)
  if (outputParameters.device == paNoDevice) {
    LA_LastError = paNoDevice;
    return LA_ERROR;
  }

  outputParameters.channelCount = 2;       /* stereo output */
  outputParameters.sampleFormat = paInt16; /* 16 bit int output */
  outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultHighOutputLatency;
  outputParameters.hostApiSpecificStreamInfo = NULL;

  PaError err;
  err = Pa_OpenStream(&LA_PlayStream,
		      NULL,
		      &outputParameters,
		      samprate,
		      paFramesPerBufferUnspecified,
		      paNoFlag,
		      cb,
		      cb_data );
  if( err != paNoError ) return LA_ERROR;

  err = Pa_StartStream(LA_PlayStream);
  if (err != paNoError) return LA_ERROR;

  return LA_OK;
}

/**
 * Play sound from a buffer of audio data.
 *
 * @param buffer A buffer of 16-bit ints
 * @param bufsize Size of the \a buffer
 * @param samprate Sampling rate of the audio
 *
 * @returns A value of type ::LaError_t where LA_OK indicates there were no problems, and a value of LA_ERROR
 * indicates that some error occurred. Use ::LA_print_last_error() to print the error.
 * 
 * @note This function is non-blocking and will return immediately (i.e. before the sound has 
 * finished playing. So use ::LA_is_playing() to determine when the sound is finished playing.
 *
 */
LaError_t LA_play(const void* const buffer, unsigned long bufsize, unsigned samprate) {

  if (LA_InitDone == 0) { /* initialize portaudio system, if needed */
    int failed = LA_init();
    if (failed) return LA_ERROR;
  }
  if (LA_is_playing() == 1) return LA_ERROR;

  LA_data.audiodata = (int16_t*) buffer;
  LA_data.len = bufsize;
  LA_data.pos_r = 0;

  PaStreamParameters outputParameters;
  outputParameters.device = Pa_GetDefaultOutputDevice();
  if (outputParameters.device == paNoDevice) {
    LA_LastError = paNoDevice;
    return LA_ERROR;
  }
  //printf( "Output device # %d.\n", outputParameters.device );
  //printf( "Output LL: %g s\n", Pa_GetDeviceInfo( outputParameters.device )->defaultLowOutputLatency );
  //printf( "Output HL: %g s\n", Pa_GetDeviceInfo( outputParameters.device )->defaultHighOutputLatency );

  outputParameters.channelCount = 2;       /* stereo output */
  outputParameters.sampleFormat = paInt16; /* 16 bit int output */
  outputParameters.suggestedLatency = Pa_GetDeviceInfo( outputParameters.device )->defaultHighOutputLatency;
  outputParameters.hostApiSpecificStreamInfo = NULL;

  PaError err;
  err = Pa_OpenStream(&LA_PlayStream,
		      NULL,
		      &outputParameters,
		      samprate,
		      paFramesPerBufferUnspecified,
		      paNoFlag,
		      _PlayCallback,
		      &LA_data );
  if( err != paNoError ) return LA_ERROR;

  err = Pa_SetStreamFinishedCallback(LA_PlayStream, _PlayStreamFinished);
  if (err != paNoError) return LA_ERROR;

  err = Pa_StartStream(LA_PlayStream);
  if (err != paNoError) return LA_ERROR;

  return LA_OK;
}

/**
 * Print last error to specified stream.
 */
void LA_print_last_error(FILE *iostream) {
  fprintf(iostream, "PortAudio Error: %s\n",Pa_GetErrorText(LA_LastError));
}

/**
 * Initialize local audio system.
 * 
 * @warning This must be called before using the local audio system.
 *
 * @returns 0 on success, or 1 if there was a problem
 */
LaError_t LA_init() {
  if (LA_InitDone) {
    int ok = LA_terminate();
    if (!ok) return LA_ERROR;
  }

  LA_LastError = Pa_Initialize();
  if (LA_LastError != paNoError) {
    LA_InitDone = 0;
    return LA_ERROR;
  }
  bzero(&LA_data, sizeof(LA_data));

  LA_InitDone = 1;
  return LA_OK;
}

/**
 * Terminate the local audio system.
 *
 * @warning It is important to call this routine when you are done with the local
 * audio system.
 *
 * @returns 0 on success, or 1 if there was a problem
 */
LaError_t LA_terminate() {
  if (LA_InitDone == 0) {   /* Can't terminate if it is not initialized */
    LA_LastError = paNotInitialized;
    return LA_ERROR;
  }
  LA_LastError = Pa_Terminate();
  if (LA_LastError != paNoError) {
    LA_InitDone = 0;
    return LA_ERROR;
  }
  LA_InitDone = 1;
  return LA_OK;
}

