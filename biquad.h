/**
 * @file biquad.h
 * @brief Biquad (second-order IIR) filter interface.
 *
 * A biquad filter is the most common building block in audio equalizers.
 * "Biquad" is short for "biquadratic" -- it describes the ratio of two
 * quadratic (degree-2) polynomials that defines the filter's transfer
 * function in the z-domain:
 *
 *   H(z) = (b0 + b1*z^-1 + b2*z^-2) / (a0 + a1*z^-1 + a2*z^-2)
 *
 * Every sample passes through five multiplications and four additions,
 * making it very efficient.  By choosing different coefficient values,
 * the same structure can produce low-pass, high-pass, band-pass,
 * notch, peaking EQ, and shelving filters.
 *
 * Based on Robert Bristow-Johnson's "Audio EQ Cookbook":
 *   https://www.w3.org/2011/audio/audio-eq-cookbook.html
 */
#ifndef BIQUAD_H
#define BIQUAD_H

#include <math.h>
#include <stdlib.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* The sample type used throughout.  Change to float if you want
 * single-precision processing (saves memory, slightly less accurate). */
typedef double smp_type;

/**
 * State and coefficients for a single biquad filter section.
 *
 * The five a-coefficients encode the filter transfer function:
 *   y[n] = a0*x[n] + a1*x[n-1] + a2*x[n-2] - a3*y[n-1] - a4*y[n-2]
 *
 * x1, x2 hold the two most recent input samples ("delay line").
 * y1, y2 hold the two most recent output samples.
 */
typedef struct {
    smp_type a0, a1, a2, a3, a4;  /* filter coefficients   */
    smp_type x1, x2;              /* input delay line       */
    smp_type y1, y2;              /* output delay line      */
} biquad;

/**
 * Process a single sample through the filter and return the result.
 * Call this once per sample in your audio loop.
 */
smp_type BiQuad(smp_type sample, biquad *b);

/**
 * Create and initialise a new biquad filter.
 *
 * @param type       One of the FILT_TYPE values (LPF, HPF, etc.).
 * @param dbGain     Gain in dB (only used for PEQ, LSH, HSH).
 * @param freq       Centre / corner frequency in Hz.
 * @param srate      Sampling rate in Hz (e.g. 44100, 48000).
 * @param bandwidth  Bandwidth in octaves (controls the filter's "width").
 * @return           A heap-allocated biquad, or nullptr on failure.
 *                   The caller must free() it when done.
 */
biquad *BiQuad_new(int type, smp_type dbGain,
                   smp_type freq, smp_type srate,
                   smp_type bandwidth);

/**
 * Filter types.
 *
 * Each type shapes the frequency response differently:
 *
 *   LPF   -- passes frequencies below the cutoff, attenuates above.
 *   HPF   -- passes frequencies above the cutoff, attenuates below.
 *   BPF   -- passes a band of frequencies, attenuates both sides.
 *   NOTCH -- removes a narrow band of frequencies (inverse of BPF).
 *   PEQ   -- boosts or cuts a band (parametric EQ).
 *   LSH   -- boosts or cuts everything below a frequency (low shelf).
 *   HSH   -- boosts or cuts everything above a frequency (high shelf).
 */
enum FILT_TYPE {
    LPF,   /* Low-pass filter   */
    HPF,   /* High-pass filter  */
    BPF,   /* Band-pass filter  */
    NOTCH, /* Notch filter      */
    PEQ,   /* Peaking band EQ   */
    LSH,   /* Low shelf filter  */
    HSH    /* High shelf filter */
};

#endif /* BIQUAD_H */
