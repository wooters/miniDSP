/**
 * @file sinedriver.c
 *
 * Example code showing how to run GCC-PHAT in order to determine the delay between several signals.
 * @ingroup GCC
 */
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include "minidsp.h"

/**
 * A simple function to rotate a vector of doubles by a given amount.
 */
void rotate(const double* const in, double* out, unsigned int n, int rot) {
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
  int    true_delays[3] = {-7,5,-2};  /* The "true" delay values. This is the answer we are looking for. */
  double amp = 10.0;       /* amplitude of the sine wave - arbitrary */
  int    sampFreq = 512;   /* frequency of the sine wave - arbitrary */

  double sampPer = 1.0/(double)sampFreq; /* sample period of the test signals */
  int    nsamps = 8192;                  /* number of samples in the test signals */
  double* siga = malloc(nsamps*sizeof(double));      /* storage for the test signal (sine wave) */
  double* sigb = malloc(nsamps*sizeof(double));      /* shifted version of test signal */
  double* sigc = malloc(nsamps*sizeof(double));      /* shifted version of test signal */
  double* sigd = malloc(nsamps*sizeof(double));      /* shifted version of test signal */
  const double** const all = malloc(4*sizeof(double*));
  all[0]=siga; all[1]=sigb; all[2]=sigc; all[3]=sigd;

  int results[3];

  /* Generate a test signal - a sine wave */
  for (unsigned i=0;i<nsamps;i++) {
    siga[i] = amp * sin((2.0 * M_PI * sampPer * i));
  }

  /* Generate delayed versions of the original sine wave */
  rotate(siga, sigb, nsamps, true_delays[0]);
  rotate(siga, sigc, nsamps, true_delays[1]);
  rotate(siga, sigd, nsamps, true_delays[2]);

  /* Compute GCC-PHAT. The GCC-PHAT values at different lags will go into lagvals[] */
  unsigned margin = 50;
  MD_get_multiple_delays(all, 4, nsamps, margin, PHAT, results); /* SIMP or PHAT */

  printf("delay A/B: %d, expected %d\n",results[0],true_delays[0]);
  printf("delay A/C: %d, expected %d\n",results[1],true_delays[1]);
  printf("delay A/D: %d, expected %d\n",results[2],true_delays[2]);

  free(all);
  free(sigd);
  free(sigc);
  free(sigb);
  free(siga);
}
