/**
 * @file minidsp_core.c
 * @brief Stateless signal math: dot product, energy, power, entropy,
 *        scaling, AGC, time-domain pitch estimation, simple effects,
 *        and window function generation.
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
 * Public API: signal analysis
 * -----------------------------------------------------------------------*/

/** Minimum normalised autocorrelation peak height accepted as voiced F0. */
#define MD_F0_ACF_PEAK_THRESHOLD 0.15

/** Three-point parabolic refinement around a discrete peak index.
 *  Returns a fractional offset in [-0.5, 0.5]. */
static double md_parabolic_offset(double y_left, double y_mid, double y_right)
{
    double denom = y_left - 2.0 * y_mid + y_right;
    if (fabs(denom) < 1e-12) return 0.0;

    double delta = 0.5 * (y_left - y_right) / denom;
    if (delta < -0.5) delta = -0.5;
    if (delta >  0.5) delta =  0.5;
    return delta;
}

/**
 * Root mean square: the standard measure of signal "loudness".
 *
 *   RMS = sqrt(1/N * sum(x[n]^2)) = sqrt(power)
 */
double MD_rms(const double *a, unsigned N)
{
    assert(a != nullptr);
    assert(N > 0);
    return sqrt(MD_power(a, N));
}

/**
 * Zero-crossing rate: fraction of adjacent sample pairs that
 * differ in sign.  Returns a value in [0.0, 1.0].
 *
 * Zero is treated as non-negative (standard convention).
 */
double MD_zero_crossing_rate(const double *a, unsigned N)
{
    assert(a != nullptr);
    assert(N > 1);
    unsigned crossings = 0;
    for (unsigned i = 1; i < N; i++) {
        if ((a[i] < 0.0) != (a[i - 1] < 0.0))
            crossings++;
    }
    return (double)crossings / (double)(N - 1);
}

/**
 * Normalised autocorrelation for lags 0..max_lag-1.
 *
 *   out[tau] = (1/R[0]) * sum_{n=0}^{N-1-tau} x[n] * x[n+tau]
 *
 * Reuses MD_dot() for the inner product.  Silent signals (R[0]=0)
 * produce all-zero output.
 */
void MD_autocorrelation(const double *a, unsigned N,
                        double *out, unsigned max_lag)
{
    assert(a != nullptr);
    assert(out != nullptr);
    assert(N > 0);
    assert(max_lag > 0 && max_lag < N);

    double r0 = MD_energy(a, N);
    if (r0 == 0.0) {
        memset(out, 0, max_lag * sizeof(double));
        return;
    }
    for (unsigned tau = 0; tau < max_lag; tau++) {
        out[tau] = MD_dot(a, a + tau, N - tau) / r0;
    }
}

/**
 * Peak detection: find local maxima above a threshold.
 *
 * A sample a[i] is a peak if:
 *   - a[i] > a[i-1]  AND  a[i] > a[i+1]  (strictly greater than neighbours)
 *   - a[i] >= threshold
 *   - distance from the last accepted peak >= min_distance
 *
 * Endpoints (i=0, i=N-1) are never peaks.
 */
void MD_peak_detect(const double *a, unsigned N, double threshold,
                    unsigned min_distance, unsigned *peaks_out,
                    unsigned *num_peaks_out)
{
    assert(a != nullptr);
    assert(peaks_out != nullptr);
    assert(num_peaks_out != nullptr);
    assert(min_distance >= 1);

    unsigned count = 0;
    unsigned last_peak = 0;

    for (unsigned i = 1; i + 1 < N; i++) {
        if (a[i] > a[i - 1] && a[i] > a[i + 1] && a[i] >= threshold) {
            if (count == 0 || i - last_peak >= min_distance) {
                peaks_out[count++] = i;
                last_peak = i;
            }
        }
    }
    *num_peaks_out = count;
}

double MD_f0_autocorrelation(const double *signal, unsigned N,
                             double sample_rate,
                             double min_freq_hz, double max_freq_hz)
{
    assert(signal != nullptr);
    assert(N >= 2);
    assert(sample_rate > 0.0);
    assert(min_freq_hz > 0.0);
    assert(max_freq_hz > min_freq_hz);

    unsigned lag_min = (unsigned)floor(sample_rate / max_freq_hz);
    unsigned lag_max = (unsigned)ceil(sample_rate / min_freq_hz);

    if (lag_min < 1) lag_min = 1;
    if (lag_max > N - 2) lag_max = N - 2;
    if (lag_min > lag_max) return 0.0;
    if (lag_max < 2) return 0.0;  /* need lag-1 and lag+1 neighbors */

    double *acf = malloc((lag_max + 1) * sizeof(double));
    if (!acf) return 0.0;

    MD_autocorrelation(signal, N, acf, lag_max + 1);

    unsigned start = (lag_min < 1) ? 1 : lag_min;
    unsigned stop  = (lag_max > N - 2) ? (N - 2) : lag_max;

    double best_peak = -DBL_MAX;
    unsigned best_lag = 0;

    for (unsigned lag = start; lag <= stop; lag++) {
        double y = acf[lag];
        if (y < MD_F0_ACF_PEAK_THRESHOLD) continue;
        if (y <= acf[lag - 1] || y <= acf[lag + 1]) continue;

        if (y > best_peak) {
            best_peak = y;
            best_lag = lag;
        }
    }

    if (best_lag == 0) {
        free(acf);
        return 0.0;
    }

    double lag_est = (double)best_lag;
    if (best_lag > 0 && best_lag + 1 <= lag_max) {
        lag_est += md_parabolic_offset(acf[best_lag - 1],
                                       acf[best_lag],
                                       acf[best_lag + 1]);
    }
    free(acf);

    if (lag_est <= 0.0) return 0.0;
    return sample_rate / lag_est;
}

/**
 * Signal mixing: weighted sum of two signals.
 *
 *   out[n] = w_a * a[n] + w_b * b[n]
 *
 * In-place safe (out may alias a or b).
 */
void MD_mix(const double *a, const double *b, double *out,
            unsigned N, double w_a, double w_b)
{
    assert(a != nullptr);
    assert(b != nullptr);
    assert(out != nullptr);
    for (unsigned i = 0; i < N; i++) {
        out[i] = w_a * a[i] + w_b * b[i];
    }
}

/* -----------------------------------------------------------------------
 * Public API: simple effects
 * -----------------------------------------------------------------------*/

void MD_delay_echo(const double *in, double *out, unsigned N,
                   unsigned delay_samples, double feedback,
                   double dry, double wet)
{
    assert(in != nullptr);
    assert(out != nullptr);
    assert(N > 0);
    assert(delay_samples > 0);
    assert(fabs(feedback) < 1.0);

    double *delay = calloc(delay_samples, sizeof(double));
    assert(delay != nullptr);

    unsigned idx = 0;
    for (unsigned n = 0; n < N; n++) {
        double x = in[n];
        double d = delay[idx];
        out[n] = dry * x + wet * d;
        delay[idx] = x + feedback * d;
        idx++;
        if (idx == delay_samples) idx = 0;
    }

    free(delay);
}

void MD_tremolo(const double *in, double *out, unsigned N,
                double rate_hz, double depth, double sample_rate)
{
    assert(in != nullptr);
    assert(out != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);
    assert(rate_hz >= 0.0);
    assert(depth >= 0.0 && depth <= 1.0);

    double phase_step = 2.0 * M_PI * rate_hz / sample_rate;
    for (unsigned n = 0; n < N; n++) {
        double lfo = 0.5 * (1.0 + sin(phase_step * (double)n));
        double gain = (1.0 - depth) + depth * lfo;
        out[n] = in[n] * gain;
    }
}

void MD_comb_reverb(const double *in, double *out, unsigned N,
                    unsigned delay_samples, double feedback,
                    double dry, double wet)
{
    assert(in != nullptr);
    assert(out != nullptr);
    assert(N > 0);
    assert(delay_samples > 0);
    assert(fabs(feedback) < 1.0);

    double *comb = calloc(delay_samples, sizeof(double));
    assert(comb != nullptr);

    unsigned idx = 0;
    for (unsigned n = 0; n < N; n++) {
        double x = in[n];
        double delayed = comb[idx];
        double y_comb = x + feedback * delayed;
        comb[idx] = y_comb;
        out[n] = dry * x + wet * y_comb;
        idx++;
        if (idx == delay_samples) idx = 0;
    }

    free(comb);
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
    assert(out != nullptr);
    assert(n > 0);

    if (n == 1) {
        out[0] = 1.0;
        return;
    }

    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 0.5 * (1.0 - cos(2.0 * M_PI * (double)i / n_minus_1));
    }
}

/**
 * Generate a Hamming window.
 *
 * The formula is:  w[i] = 0.54 - 0.46 * cos(2*pi*i / (n-1))
 */
void MD_Gen_Hamming_Win(double *out, unsigned n)
{
    assert(out != nullptr);
    assert(n > 0);

    if (n == 1) {
        out[0] = 1.0;
        return;
    }

    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 0.54 - 0.46 * cos(2.0 * M_PI * (double)i / n_minus_1);
    }
}

/**
 * Generate a Blackman window.
 *
 * The formula is:
 *   w[i] = 0.42
 *        - 0.5  * cos(2*pi*i / (n-1))
 *        + 0.08 * cos(4*pi*i / (n-1))
 */
void MD_Gen_Blackman_Win(double *out, unsigned n)
{
    assert(out != nullptr);
    assert(n > 0);

    if (n == 1) {
        out[0] = 1.0;
        return;
    }

    double n_minus_1 = (double)(n - 1);
    for (unsigned i = 0; i < n; i++) {
        double p = 2.0 * M_PI * (double)i / n_minus_1;
        out[i] = 0.42 - 0.5 * cos(p) + 0.08 * cos(2.0 * p);
    }
}

/**
 * Generate a rectangular window (all ones).
 */
void MD_Gen_Rect_Win(double *out, unsigned n)
{
    assert(out != nullptr);
    assert(n > 0);
    for (unsigned i = 0; i < n; i++) {
        out[i] = 1.0;
    }
}
