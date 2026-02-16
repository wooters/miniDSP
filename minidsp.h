/**
 * @file minidsp.h
 * @brief A mini library of DSP (Digital Signal Processing) routines.
 *
 * This header declares functions for:
 *   - Basic signal measurements (energy, power, entropy)
 *   - Signal scaling and gain adjustment
 *   - Window generation (Hanning window)
 *   - Generalized Cross-Correlation (GCC-PHAT) for delay estimation
 *
 * These are the kinds of building blocks you'd use in an audio processing
 * pipeline -- for example, estimating which direction a sound came from
 * using a pair of microphones.
 */
#ifndef MINIDSP_H
#define MINIDSP_H

#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/* M_PI is not guaranteed by C11.  Define it if the system didn't. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#include <float.h>
#include <complex.h>
#include <fftw3.h>

/* -----------------------------------------------------------------------
 * Basic signal measurement functions
 * -----------------------------------------------------------------------*/

/**
 * Compute the dot product of two vectors.
 * The dot product is the sum of element-wise products: a[0]*b[0] + a[1]*b[1] + ...
 */
double MD_dot(const double *a, const double *b, unsigned N);

/**
 * Compute the normalized entropy of a distribution.
 * Returns a value between 0.0 (all energy concentrated in one bin)
 * and 1.0 (energy spread equally across all bins).
 *
 * @param clip  If true, ignore negative values. If false, square all values first.
 */
double MD_entropy(const double *a, unsigned N, bool clip);

/** Compute signal energy: sum of squared samples. */
double MD_energy(const double *a, unsigned N);

/** Compute signal power: energy divided by the number of samples. */
double MD_power(const double *a, unsigned N);

/** Compute signal power in decibels: 10 * log10(power). */
double MD_power_db(const double *a, unsigned N);

/* -----------------------------------------------------------------------
 * Signal scaling and conditioning
 * -----------------------------------------------------------------------*/

/**
 * Map a single value from one range to another.
 * Example: MD_scale(5, 0, 10, 0, 100) returns 50.
 */
double MD_scale(double in,
                double oldmin, double oldmax,
                double newmin, double newmax);

/** Map every element of a vector from one range to another. */
void MD_scale_vec(double *in, double *out, unsigned N,
                  double oldmin, double oldmax,
                  double newmin, double newmax);

/**
 * Fit values within [newmin, newmax].
 * If all values already fit, they are copied unchanged.
 * Otherwise the entire vector is rescaled.
 */
void MD_fit_within_range(double *in, double *out, unsigned N,
                         double newmin, double newmax);

/**
 * Automatic Gain Control: scale a signal so its power matches
 * the requested dB level, then clip to [-1, 1].
 */
void MD_adjust_dblevel(const double *in, double *out,
                       unsigned N, double dblevel);

/* -----------------------------------------------------------------------
 * Window generation
 * -----------------------------------------------------------------------*/

/**
 * Generate a Hanning window of length n.
 * A Hanning window tapers the edges of a signal to zero, which reduces
 * spectral leakage when you later take an FFT.
 */
void MD_Gen_Hann_Win(double *out, unsigned n);

/* -----------------------------------------------------------------------
 * Resource cleanup
 * -----------------------------------------------------------------------*/

/** Free all internally cached FFT plans and buffers. Call when done. */
void MD_shutdown(void);

/* -----------------------------------------------------------------------
 * Generalized Cross-Correlation (GCC) for delay estimation
 * -----------------------------------------------------------------------
 *
 * Given two microphone signals that captured the same sound source,
 * GCC-PHAT estimates how many samples one signal is delayed relative
 * to the other.  This is the basis of acoustic source localisation.
 *
 * The algorithm:
 *   1. FFT both signals.
 *   2. Multiply one spectrum by the conjugate of the other (cross-spectrum).
 *   3. Apply a weighting (PHAT normalises by magnitude, sharpening the peak).
 *   4. Inverse-FFT back to the time domain.
 *   5. The position of the peak tells you the delay in samples.
 */

/** Weighting types for Generalized Cross-Correlation. */
enum MD_GCC_WEIGHTING_TYPE {
    SIMP, /**< Simple 1/N weighting (basic cross-correlation) */
    PHAT  /**< Phase Transform weighting (sharper peaks, more robust to noise) */
};

/**
 * Estimate delays between a reference signal and M-1 other signals.
 *
 * @param sigs        Array of M pointers to signals (sigs[0] is the reference).
 * @param M           Number of signals.
 * @param N           Length of each signal (all must be the same length).
 * @param margin      Search +/- this many samples around zero-lag.
 * @param weightfunc  SIMP or PHAT (see ::MD_GCC_WEIGHTING_TYPE).
 * @param outdelays   Output array of M-1 delay values (must be pre-allocated).
 */
void MD_get_multiple_delays(const double **sigs, unsigned M, unsigned N,
                            unsigned margin, int weightfunc,
                            int *outdelays);

/**
 * Estimate the delay between two signals.
 *
 * @param siga        First signal.
 * @param sigb        Second signal.
 * @param N           Length of both signals.
 * @param ent         If non-NULL, receives the normalised entropy of the
 *                    correlation peak region (closer to 1.0 = less trustworthy).
 * @param margin      Search +/- this many samples around zero-lag.
 * @param weightfunc  SIMP or PHAT.
 * @return            Delay in samples (positive = sigb lags siga).
 */
int MD_get_delay(const double *siga, const double *sigb, unsigned N,
                 double *ent, unsigned margin, int weightfunc);

/**
 * Compute the full generalized cross-correlation between two signals.
 *
 * @param siga        First signal.
 * @param sigb        Second signal.
 * @param N           Length of both signals.
 * @param lagvals     Output array of N doubles (pre-allocated).
 *                    The zero-lag value is at index ceil(N/2).
 * @param weightfunc  SIMP or PHAT.
 */
void MD_gcc(const double *siga, const double *sigb, unsigned N,
            double *lagvals, int weightfunc);

#endif /* MINIDSP_H */
