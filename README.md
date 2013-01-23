#MiniDSP

A small C-lib containing a few basic DSP routines, and other helpful functions for dealing with audio.

##Dependencies

* FFTW - http://www.fftw.org/ (TODO: allow for use of other FFT libs)
* Portaudio - http://portaudio.com/ (for live audio input/output)
* libsndfile - http://www.mega-nerd.com/libsndfile/ (for file-based input/output)
* gnuplot_i - http://ndevilla.free.fr/gnuplot/ (for plotting with gnuplot, assumes you have gnuplot installed)

##Examples
### Record from the mic at a rate of 16kHz and plot the data with gnuplot:

```c
/**
 * @file testliverec.c
 */ 
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "gnuplot_i.h"
#include "liveio.h"

int main(void) {
  unsigned sampRate = 16000;
  unsigned nsamps = sampRate * 5;
  int16_t* buffer = (int16_t*) calloc(nsamps, sizeof(int16_t));

  /* Record into the buffer */  
  LA_record(buffer, nsamps, sampRate, LA_REC_ONCE); 
  printf("Recording...\n");fflush(stdout);
  while(LA_is_recording() == 1) sleep(1);

  /* Terminate the audio recording session */
  LA_terminate();

  /* Convert 16-bit ints to doubles for plotting */
  double* dbuf = calloc(nsamps, sizeof(double));
  for (unsigned i=0;i<nsamps;i++) dbuf[i] = (double) buffer[i];

  /* Plot with gnuplot */
  gnuplot_ctrl* h = gnuplot_init();
  gnuplot_setstyle(h,"lines");
  gnuplot_plot_x(h,dbuf,nsamps,"Live Audio");
  gnuplot_close(h);

  /* Clean up */
  free(dbuf);
  free(buffer);

  return 0;
}
```

Compile with:

```make
LMDSP= $(HOME)/src/miniDSP/libminidsp.a
GP= $(HOME)/src/gnuplot_i/src
GP_OBJ= $(GP)/gnuplot_i.o

CFLAGS = -O3 -Wall

%.o:%.c
	$(CC) $(CFLAGS) -I $(GP) -I .. -c $*.c -o $@ 

testliverec: testliverec.o $(GP_OBJ) $(LMDSP)
	$(CC) $(CFLAGS) -o $@ $(TLR_OBJS) $(LDFLAGS) -L$(HOME)/src/miniDSP -lminidsp -lportaudio

```


##License

This open-source software is licensed under the MIT License. See the
LICENSE file for details.

## Author
Chuck Wooters - <wooters@icsi.berkeley.edu>
