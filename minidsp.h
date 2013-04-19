/**
 * @file minidsp.h
 * @brief A mini library of DSP-related routines.
 *
 */
#ifndef MINIDSP_H
#define MINIDSP_H

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <complex.h>
#include <fftw3.h>

double MD_dot(const double* const a, const double* const b, const unsigned N);

double MD_entropy(const double* const a, const unsigned N, const bool clip);
double MD_energy(const double* const a, const unsigned N);
double MD_power(const double* const a, const unsigned N);
double MD_power_db(const double* const a, const unsigned N);


double MD_scale(const double in, const double oldmin, const double oldmax,
		const double newmin, const double newmax);
void   MD_scale_vec(double* const in, double* const out, const unsigned N, 
		    const double oldmin, const double oldmax,
		    const double newmin, const double newmax);
void   MD_fit_within_range(double* const in, double* const out, 
			   const unsigned N, const double newmin, 
			   const double newmax);
void   MD_adjust_dblevel(const double* const in, double* const out, 
			 const unsigned N, const double dblevel);

void   MD_Gen_Hann_Win(double* out, unsigned n);

void   MD_shutdown(void);

/** 
 * Possible weighting types for Generalized Cross Correlation.
 *
 */
enum MD_GCC_WEIGHTING_TYPE { 
  SIMP, /**< Simple \f$1/N\f$ weighting */
  PHAT  /**< Phase Transform (PHAT) weighting */
};

/* 
 * Routines for computing delays between signals 
 */
void MD_get_multiple_delays(const double** const sigs, const unsigned M, const unsigned N, 
			    const unsigned margin, const int weightfunc, int* outdelays);
 
int MD_get_delay(const double* const siga, const double* const sigb, const unsigned N,
		 double* const ent,
		 const unsigned margin, 
		 const int weightfunc);

void MD_gcc(const double* const siga, const double* const sigb, const unsigned N,
	    double* const lagvals, const int weightfunc);


#endif


