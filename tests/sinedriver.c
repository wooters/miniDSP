/**
 * @file sinedriver.c
 *
 * Example code showing how to run GCC-PHAT in order to determine the delay between two signals.
 * @ingroup GCC
 */
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "minidsp.h"
#include "gnuplot_i.h"

/**
 * A simple function to rotate a vector of doubles by a given amount.
 */
void rotate(double* in, double* out, unsigned int n, int rot) {
  for (unsigned int i=0; i<n ;i++) {
    unsigned int j = (i+rot)%n;
    out[j] = in[i];
  }
}

/**
 * A simple program to test the GCC-PHAT routine.
 *
 * @see ::get_delay()
 */
int main()
{
  int    true_delay = -7;  /* The "true" delay value. This is the answer we are looking for. */
  double amp = 10.0;       /* amplitude of the sine wave - arbitrary */
  int    sampFreq = 512;   /* frequency of the sine wave - arbitrary */

  double sampPer = 1.0/(double)sampFreq; /* sample period of the test signals */
  int    nsamps = 8192;                  /* number of samples in the test signals */
  double siga[nsamps];                   /* storage for the test signal (sine wave) */
  double sigb[nsamps];                   /* shifted version of test signal */

  /* Generate a test signal - a sine wave */
  for (unsigned i=0;i<nsamps;i++) {
    siga[i] = amp * sin((2.0 * M_PI * sampPer * i));
  }

  /* Generate a delayed version of the sine wave */
  rotate(siga, sigb, nsamps, true_delay);

  /* Compute GCC-PHAT. The GCC-PHAT values at different lags will go into lagvals[] */
  unsigned margin = 50;
  int delay = MD_get_delay(siga, sigb, nsamps, margin, PHAT); /* SIMP or PHAT */

  printf("gcc:\n");
  printf("  estimated delay: %d\n",delay);

  /* Set up for plotting with GNUplot */
  //gnuplot_ctrl* h = gnuplot_init();
  //gnuplot_setstyle(h,"lines");
  /*gnuplot_cmd(h,"set arrow from %d,%d to %d,%d",50,(int)maxval/2,maxindex,(int)maxval);*/
  //gnuplot_plot_x(h,lagvals,nsamps,"GCC-PHAT Cross Correlation");
  //gnuplot_close(h); 
}
