/**
 * @file minidsp_core.c
 * @brief Stateless signal math: dot product, energy, power, entropy,
 *        scaling, AGC, and Hanning window generation.
 * @author Chuck Wooters <wooters@hey.com>
 * @copyright 2013 International Computer Science Institute
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Public API: basic signal measurements
 * -----------------------------------------------------------------------*/

/**
 * Dot product of two vectors a and b, each of length N.
 *
 *   dot = a[0]*b[0] + a[1]*b[1] + ... + a[N-1]*b[N-1]
 *
 * This is one of the most fundamental operations in DSP.  The dot
 * product measures how "similar" two signals are -- if they are
 * identical the result is large; if they are unrelated it is near zero.
 */
double MD_dot(const double *a, const double *b, unsigned N)
{
    double d = 0.0;
    for (unsigned i = 0; i < N; i++) {
        d += a[i] * b[i];
    }
    return d;
}

/**
 * Signal energy: the sum of squared samples.
 *
 *          N-1
 *   E  =  SUM  a[n]^2
 *          n=0
 *
 * Energy tells you "how loud" a signal is overall.  A silent signal
 * has zero energy; a loud one has high energy.  Notice that squaring
 * makes all values positive, so negative samples contribute too.
 */
double MD_energy(const double *a, unsigned N)
{
    assert(a != nullptr);
    if (N == 1) return a[0] * a[0];
    return MD_dot(a, a, N);
}

/**
 * Signal power: energy divided by the number of samples.
 *
 *   P = E / N
 *
 * Power is the "average energy per sample".  It lets you compare
 * signals of different lengths on an equal footing.
 */
double MD_power(const double *a, unsigned N)
{
    assert(a != nullptr);
    assert(N > 0);
    return MD_energy(a, N) / (double)N;
}

/**
 * Signal power expressed in decibels (dB).
 *
 *   P_dB = 10 * log10(P)
 *
 * Decibels are a logarithmic scale commonly used in audio engineering.
 * Every +10 dB means the power is 10x larger.  A floor of 1e-10 is
 * used to avoid log(0), which would be negative infinity.
 */
double MD_power_db(const double *a, unsigned N)
{
    assert(a != nullptr);
    assert(N > 0);
    double p = fmax(1.0e-10, MD_power(a, N));
    return 10.0 * log10(p);
}

/* -----------------------------------------------------------------------
 * Public API: signal scaling and conditioning
 * -----------------------------------------------------------------------*/

/**
 * Linearly map a value from one range to another.
 *
 * Formula:  out = (in - oldmin) * (newmax - newmin) / (oldmax - oldmin) + newmin
 *
 * For example, mapping a value of 5 from the range [0, 10] into [0, 100]
 * gives 50.  This is sometimes called "lerp" (linear interpolation).
 */
double MD_scale(double in,
                double oldmin, double oldmax,
                double newmin, double newmax)
{
    return (in - oldmin) * (newmax - newmin) / (oldmax - oldmin) + newmin;
}

/**
 * Apply MD_scale() to every element of a vector.
 */
void MD_scale_vec(double *in, double *out, unsigned N,
                  double oldmin, double oldmax,
                  double newmin, double newmax)
{
    assert(in != nullptr);
    assert(out != nullptr);
    assert(oldmin < oldmax);
    assert(newmin < newmax);
    if (N == 0) return;

    double scale = (newmax - newmin) / (oldmax - oldmin);
    for (unsigned i = 0; i < N; i++) {
        out[i] = (in[i] - oldmin) * scale + newmin;
    }
}

/**
 * Squeeze values into [newmin, newmax] only if they don't already fit.
 *
 * If the input range already lies inside [newmin, newmax], the values
 * are copied as-is (no stretching).  Otherwise, the full input range
 * is mapped onto the new range.
 */
void MD_fit_within_range(double *in, double *out, unsigned N,
                         double newmin, double newmax)
{
    assert(in != nullptr);
    assert(out != nullptr);
    assert(newmin < newmax);
    if (N == 0) return;

    /* Find the actual min and max of the input */
    double in_min = in[0];
    double in_max = in[0];
    for (unsigned i = 1; i < N; i++) {
        if (in[i] < in_min) in_min = in[i];
        if (in[i] > in_max) in_max = in[i];
    }

    if (in_min >= newmin && in_max <= newmax) {
        /* Already fits -- just copy without scaling */
        for (unsigned i = 0; i < N; i++) {
            out[i] = in[i];
        }
    } else {
        MD_scale_vec(in, out, N, in_min, in_max, newmin, newmax);
    }
}

/**
 * Automatic Gain Control (AGC): adjust a signal to a target dB level.
 *
 * The input signal is assumed to be in the range [-1.0, 1.0].
 * The function computes a gain factor so that the output signal
 * has the desired power level (in dB), then clips to [-1.0, 1.0]
 * if any sample exceeds that range.
 *
 * Based on Jaydeep Dhole's AGC implementation for MATLAB.
 */
void MD_adjust_dblevel(const double *in, double *out,
                       unsigned N, double dblevel)
{
    /* Convert the target dB level back to linear power */
    double desired_power = pow(10.0, dblevel / 10.0);

    /* Compute the gain factor:
     * We want:  (1/N) * sum(out[i]^2) = desired_power
     * Since out[i] = in[i] * gain:
     *   (gain^2 / N) * sum(in[i]^2) = desired_power
     *   gain = sqrt(desired_power * N / energy_in) */
    double input_energy = MD_energy(in, N);
    double gain = sqrt((desired_power * (double)N) / input_energy);

    /* Apply the gain and check for out-of-range values */
    bool out_of_range = false;
    for (unsigned i = 0; i < N; i++) {
        out[i] = in[i] * gain;
        if (out[i] > 1.0 || out[i] < -1.0) {
            out_of_range = true;
        }
    }

    /* If any samples exceed [-1.0, 1.0], rescale everything to fit */
    if (out_of_range) {
        MD_fit_within_range(out, out, N, -1.0, 1.0);
    }
}

/**
 * Compute the normalised entropy of a distribution.
 *
 * Entropy measures how "spread out" a distribution is:
 *   - 0.0 means all the energy is in a single bin (very peaky).
 *   - 1.0 means the energy is spread equally (completely flat).
 *
 * The Shannon entropy formula is:  H = -SUM( p_i * log2(p_i) )
 * We normalise by dividing by log2(N) so the result is always in [0, 1].
 *
 * @param a     Array of values representing the distribution.
 * @param N     Length of the array.
 * @param clip  If true, ignore negative values (treat them as 0).
 *              If false, use a[i]^2 so all values become non-negative.
 * @return      Normalised entropy in [0.0, 1.0].
 */
double MD_entropy(const double *a, unsigned N, bool clip)
{
    assert(a != nullptr);

    if (N <= 1) return 0.0;

    /* Maximum possible entropy for N bins (uniform distribution) */
    double max_entropy = log2((double)N);

    /* First pass: compute the total (for normalisation into probabilities) */
    double total = 0.0;
    if (clip) {
        for (unsigned i = 0; i < N; i++) {
            total += (a[i] < 0.0) ? 0.0 : a[i];
        }
    } else {
        for (unsigned i = 0; i < N; i++) {
            total += a[i] * a[i];
        }
    }

    /* If the total is zero (e.g., all-zeros input), entropy is undefined.
     * We return 0.0 because there is no meaningful distribution. */
    if (total == 0.0) return 0.0;

    /* Second pass: compute H = -SUM( p_i * log2(p_i) ) */
    double entropy = 0.0;
    for (unsigned i = 0; i < N; i++) {
        if (a[i] == 0.0) continue;
        if (clip && a[i] < 0.0) continue;

        double p;
        if (clip) {
            p = a[i] / total;
        } else {
            p = (a[i] * a[i]) / total;
        }

        entropy += p * log2(p);  /* Note: p*log2(p) is always <= 0 */
    }

    /* Negate and normalise so the result is in [0, 1] */
    return -entropy / max_entropy;
}

/* -----------------------------------------------------------------------
 * Public API: window generation
 * -----------------------------------------------------------------------*/

/**
 * Generate a Hanning (raised cosine) window.
 *
 * The formula is:  w[i] = 0.5 * (1 - cos(2*pi*i / (n-1)))
 *
 * Windowing is essential before taking an FFT because real-world
 * signals don't start and end at exactly zero.  Without a window,
 * the abrupt edges create false high-frequency content ("spectral
 * leakage").  The Hanning window smoothly tapers the signal to zero
 * at both ends, reducing this artefact.
 */
void MD_Gen_Hann_Win(double *out, unsigned n)
{
    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / n_minus_1));
    }
}
