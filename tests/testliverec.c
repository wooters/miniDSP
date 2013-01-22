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
