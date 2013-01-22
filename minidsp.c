/**
 * @file minidsp.c
 * @brief A small set of DSP routines.
 * @author Chuck Wooters <wooters@icsi.berkeley.edu>
 * @copyright 2013 International Computer Science Institute
 */
/**
 * @defgroup GCC Generalized Cross-Correlation
 *
 * Assuming we have two (or more) spatially separated
 * microphones \f$r_1(t)\f$ and \f$r_2(t)\f$ and a sound source
 * \f$s(t)\f$, we can model the time delay between the two received
 * signals as \f$r_1(t) = s(t) + n_1(t)\f$ and \f$r_2(t) = s(t-D) + * n_2(t)\f$ 
 * with \f$0\le t \le T\f$ where \f$T\f$ is the observation interval,
 * \f$n_1(t)\f$ and \f$n_2(t)\f$ represent
 * additive noise at each microphone and \f$D\f$ is the time delay
 * between the two signals.
 *
 * The estimated time delay \f$D\f$ between the two signals is
 * obtained by finding the time lag that maximizes the
 * cross-correlation between the signals received at the two
 * microphones. The approach used here to find the time lag is the
 * generalized cross correlation with a phase transform weighting
 * (GCC-PHAT). The phase transform weighting helps to enhance the
 * peaks in the cross-correlation function making it easier to locate
 * the maximum.
 *
 */
#include "minidsp.h"

static unsigned      _N = 0;        /**< size of input vector */
static double* siga_loc = NULL;     /**< local copy of input A */
static double* sigb_loc = NULL;     /**< local copy of input B */
static double* lags_loc = NULL;     /**< local version of computed lag values */
static fftw_complex* ffta  = NULL;  /**< fft of input A */
static fftw_complex* fftb  = NULL;  /**< fft of input B */
static fftw_complex* xspec = NULL;  /**< the cross-spectrum of A & B */
static double* xcorr = NULL;        /**< the cross-correlation of A & B */

static fftw_plan pa = NULL; /**< FFTW plan for FFT(A) */
static fftw_plan pb = NULL; /**< FFTW plan for FFT(B) */
static fftw_plan px = NULL; /**< FFTW plan for IFFT(xspec) */

/**
 * Free up memory allocated for doing GCC
 * @ingroup GCC
 */
void _xcorr_free() {
  if (xspec)  fftw_free(xspec);
  if (fftb)   fftw_free(fftb);
  if (ffta)   fftw_free(ffta);

  if (lags_loc) free(lags_loc);
  if (xcorr)    free(xcorr);
  if (sigb_loc) free(sigb_loc);
  if (siga_loc) free(siga_loc);
}

/**
 * Allocate memory for running GCC functions
 * @ingroup GCC
 */
void _xcorr_malloc() {
  siga_loc     = calloc(_N, sizeof(double));
  sigb_loc     = calloc(_N, sizeof(double));
  xcorr        = calloc(_N+1, sizeof(double));
  lags_loc     = calloc(_N, sizeof(double));

  ffta         = fftw_alloc_complex(_N);
  fftb         = fftw_alloc_complex(_N);
  xspec        = fftw_alloc_complex(_N);
}

/**
 * Clean-up when finished running GCC.
 * @ingroup GCC
 */
void _xcorr_teardown() {
  if (pa) fftw_destroy_plan(pa);
  if (pb) fftw_destroy_plan(pb);
  if (px) fftw_destroy_plan(px);
  _xcorr_free();
  fftw_cleanup();
}

/**
 * Prepare for running GCC.
 * @ingroup GCC
 */
void _xcorr_setup() {
  _xcorr_teardown();
  _xcorr_malloc();

  pa = fftw_plan_dft_r2c_1d(_N, siga_loc, ffta,  FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  pb = fftw_plan_dft_r2c_1d(_N, sigb_loc, fftb,  FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
  px = fftw_plan_dft_c2r_1d(_N, xspec,    xcorr, FFTW_ESTIMATE|FFTW_DESTROY_INPUT);
}

/**
 * Determine if we need to run \ref ::_xcorr_setup()
 * @ingroup GCC
 */
void _gcc_setup(const unsigned N) {
  if (_N != N) {/* do we need a new vector size? */
    /* update static variables */
    _N = N;               
    _xcorr_setup();        /* set up the static vectors and FFT params */
  }
}

/**
 * Shift the output of an FFT.
 *
 * The index of the mid-point in the output will be located at: ceil(_N/2)
 * @ingroup GCC
 */
void _fftshift(const double* const in, double* const out, const unsigned N) {
  /* mid-point of out[] will be located at index ceil(N/2) */
  unsigned xx = (unsigned) floor((double) N/2.0);

  /* Copy last half of in[] to first half of out[] */
  memcpy(out,&in[xx],sizeof(double)*(N-xx));

  /* Copy first half of in[] to end of out[] */
  memcpy(&out[N-xx],in,sizeof(double)*xx);
}

/**
 * Find the index and max value in the given array of doubles.
 * @ingroup GCC
 */
void _max_index(const double* const a, const unsigned N,
	      double* const max, unsigned* const maxi) {

  if (N <= 1) {
    *max = a[0]; *maxi = 0;
  } else {
    unsigned maxi_t = 0;
    double max_t = a[0];
    for (unsigned i=1;i<N;i++){
      if (a[i] > max_t) {
	max_t = a[i]; maxi_t = i;
      }
    }
    *max = max_t;
    *maxi = maxi_t;
  }
}

/**
 * Compute the delays between a reference signal and several other
 * signals. This function performs a generalized cross-correlation
 * (via the ::gcc() function) across multiple input vectors using the
 * specified weighting function. The delays (in samples) between the
 * signals will be returned in the \a outdelays array.
 * 
 * @param sigs An array of pointers to vectors of doubles containing
 * the signals. The first signal in this array will be considered the
 * reference signal.
 * @param M The number of vectors in \a sigs
 * @param N The length of each of the vectors in \a sigs
 * @param margin Search for the delay over a window of \f$\pm\f$\a
 * margin samples centered around the 0-lag
 * @param weightfunc A ::GCC_WEIGHTING_TYPE value to specify the
 * weighting function type to use in the GCC
 * @param outdelays A pointer to an array of length \a N-1 ints in which to store the delays
 *
 * @see ::MD_gcc() ::MD_get_delay()
 * @ingroup GCC
 *
 * @warning This function assumes that the memory for \a outdelays has already
 * been allocated and that it has enough space to store \a N-1 int values.
 *
 */
void MD_get_multiple_delays(const double** const sigs, const unsigned M, const unsigned N, 
			    const unsigned margin, const int weightfunc, int* const outdelays) 
{
  if (M < 2) return;
  for (unsigned i=0;i<M-1;i++) {
    outdelays[i] = MD_get_delay(sigs[0],sigs[i+1], N, margin, weightfunc);
  }
}

/**
 * Compute the delay between a pair of signals. This function performs
 * a generalized cross-correlation (via the ::gcc() function) between
 * two input vectors using the specified weighting funciton. The delay
 * (in samples) between the two signals will be returned.
 * 
 * @param siga first vector of doubles
 * @param sigb second vector of doubles
 * @param N length of vectors \a siga and \a sigb
 * @param margin the max delay will be searched for over a window of \f$\pm\f$\a margin samples centered around the 0-lag
 * @param weightfunc a ::GCC_WEIGHTING_TYPE value to specify the weighting function type to use in the GCC
 * @return An integer value indicating the delay between \a siga and \a sigb.
 *
 * @see ::MD_gcc()
 * @ingroup GCC
 *
 */
int MD_get_delay(const double* const siga, const double* const sigb, const unsigned N,
		 const unsigned margin, const int weightfunc)
{
  _gcc_setup(N);

  /* Compute the Generalized Cross-Correlation using the given weighting function */
  MD_gcc(siga, sigb, N, lags_loc, weightfunc);

  /* 
   * Search within +- margin samples to find the best delay value 
   */

  /* First, make sure the margin is within the bounds of the computed lags */
  unsigned center_i = ceil(N/2.0); /* index of the 0-lag in lags_loc[] */
  unsigned newmargin=margin;
  if (((int)(center_i - newmargin)) < 0) {
    newmargin = center_i;
  }
  if ((center_i + newmargin) >= N) {
    newmargin = (N-1) - center_i;
  }

  /* Compute the begin index and length of the lags_loc[] array */
  unsigned start_i = center_i-newmargin;
  unsigned len = 2*newmargin+1;
  double maxval;
  unsigned max_i;

  /* Search within the desired region (+-margin) for the max lag */
  _max_index(lags_loc+start_i,len,&maxval,&max_i);

  return (int)(max_i - newmargin);
}

/**
 * Perform generalized cross-correlation. This function performs a
 * generalized cross-correlation (GCC) between two input vectors using
 * a specified weighting funciton and fills an output array with 
 * cross-correlation values. 
 * 
 * @param siga first vector of doubles
 * @param sigb second vector of doubles
 * @param N length of vectors \a siga and \a sigb 
 * @param lagvals pointer to an \a N length vector of doubles in which to store the results. The lag values will be organized such that the 0-lag is located at index ceil(\a N / 2)
 * @param weightfunc a ::GCC_WEIGHTING_TYPE value to specify the weighting function type to use in the GCC
 *
 * @ingroup GCC
 * @see ::MD_get_delay()
 *
 * @warning This function assumes that the memory for \a lagvals has been already
 * been allocated and that it has enough space to store \a N doubles.
 */

void MD_gcc(const double* const siga, const double* const sigb, const unsigned N,
	      double* const lagvals, const int weightfunc)
{
  _gcc_setup(N);

  /* copy the input arrays into the local arrays */
  memcpy(siga_loc,siga,_N*sizeof(double));
  memcpy(sigb_loc,sigb,_N*sizeof(double));

  /* Perform the FFTs on the two input signals.
   *   The output will be stored in ffta[] and fftb[]. These
   *   are arrays of complex numbers and will contain only 
   *   the non-negative frequencies, plus one element.
   */
  fftw_execute(pa);
  fftw_execute(pb);

  /* Compute GCC with the specified weighting function
   *
   * Note that we only need to loop over the first _N/2+1 items
   * in ffta[] and fft[b] since the remainder of the items
   * are all 0's (FFTW only computes the non-negative frequencies
   * for real-input FFTs.
   */
  unsigned xspec_bound = _N/2+1;
  switch (weightfunc) {
  case PHAT:
    for (unsigned i = 0; i <= xspec_bound; i++) {
      double complex tmp = fftb[i] * conj(ffta[i]); /* cross-spectra */
      xspec[i] = tmp / (cabs(tmp)+DBL_MIN); /* adding DBL_MIN to prevent div by 0 */
    }
    break;
  default: /* weightfunc != PHAT */
    for (unsigned i = 0; i <= xspec_bound; i++) {
      double complex tmp = fftb[i] * conj(ffta[i]); /* cross-spectra */
      xspec[i] = tmp / (double) _N;
    }
  }

  /* Inverse Fourier transform of the cross-spectra. Output will be placed in xcorr[] */
  fftw_execute(px);

  /*
   * Shift the values in xcorr[] so that the 0th lag is at the center of
   * the output array. 
   * [Note: the index of the center value in the output will be: ceil(_N/2) ]
   */
  _fftshift(xcorr, lagvals, _N);

}

