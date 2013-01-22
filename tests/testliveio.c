/**
 * @file testliveio.c
 *
 * This program will test the liveio functions by starting a
 * continuously-running recording from the local microphone.  When a
 * key is hit, the recording will stop and whatever is in the buffer
 * will be played out the speaker.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <inttypes.h>
#include "liveio.h"

const unsigned sampFreq = 16000;
const unsigned nSecs    = 5; /* number of seconds in buffer */

void genSine(int16_t* buf, unsigned n, unsigned freq, unsigned sampfreq) {
  for(unsigned i=0;i<n;i++) {
    buf[i] = (int16_t) floor(30000.0*sin(freq * (2.0*M_PI) * i / sampfreq));
  }
}

int main(void) {
  int failed = 0;
  unsigned int bufSize = sampFreq * nSecs;
  int16_t* buffer;

  buffer = (int16_t*) calloc(bufSize, sizeof(int16_t));

  /* record continuously into buffer[]. If it reaches the end
   * of the buffer, it will just wrap around and continue
   * recording at the beginning of the buffer.
   */
  failed = LA_record(buffer, bufSize, sampFreq, LA_REC_CONT); 
  if (failed) {
    printf("Failed to record audio.\n");
    LA_print_last_error(stdout);
    exit(-1);
  }
  printf("Recording... (hit any key to stop)\n");fflush(stdout);
  getchar();
  LA_stop_recording();


  /* Recording is done, so play whatever is in the buffer[]
   * starting at the beginning.
   */
  failed = LA_play(buffer,bufSize,sampFreq);
  if (failed) {
    printf("Failed to play audio.\n");
    LA_print_last_error(stdout);
    exit(-1);
  }
  printf("Playing...\n");fflush(stdout);
  while(LA_is_playing()==1) {
    sleep(1);
  }

  /* Terminate the local audio session */
  failed = LA_terminate();
  if (failed) {
    printf("Failed to terminate audio session.\n");
    LA_print_last_error(stdout);
    exit (failed);
  }

  free(buffer);

  exit (0);
}

