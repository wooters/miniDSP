/**
 * @file minidsp.c
 * @brief Core DSP routines: signal measurements, scaling, and GCC-PHAT.
 * @author Chuck Wooters <wooters@hey.com>
 * @copyright 2013 International Computer Science Institute
 *
 * This file implements the core digital signal processing functions
 * declared in minidsp.h.  The most interesting algorithm here is
 * GCC-PHAT (Generalized Cross-Correlation with Phase Transform),
 * which figures out the time delay between two microphone signals.
 *
 * How GCC-PHAT works (the big picture):
 *
 *   Imagine you clap your hands in a room with two microphones.
 *   Microphone A hears the clap slightly before microphone B because
 *   it is closer to you.  GCC-PHAT measures that tiny time difference
 *   by looking at how the two signals line up in the frequency domain.
 *
 *   Mathematically:
 *     1. Take the FFT of both signals:  FFT(A), FFT(B)
 *     2. Compute the cross-spectrum:    X = FFT(B) * conj(FFT(A))
 *     3. Normalise by magnitude:        X_phat = X / |X|
 *        (This "phase transform" sharpens the correlation peak.)
 *     4. Inverse-FFT back to time:      xcorr = IFFT(X_phat)
 *     5. The index of the peak in xcorr tells you the delay.
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Static (file-scope) variables for FFT caching
 *
 * Creating FFTW plans is expensive, so we cache them.  As long as the
 * signal length N stays the same between calls, we reuse the existing
 * plans and buffers.  If N changes, we tear everything down and rebuild.
 * -----------------------------------------------------------------------*/

static unsigned      _N        = 0;    /* Current cached signal length        */
static double       *siga_loc  = NULL; /* Local copy of input signal A        */
static double       *sigb_loc  = NULL; /* Local copy of input signal B        */
static double       *lags_loc  = NULL; /* Shifted cross-correlation output    */
static fftw_complex *ffta      = NULL; /* FFT of signal A                     */
static fftw_complex *fftb      = NULL; /* FFT of signal B                     */
static fftw_complex *xspec     = NULL; /* Cross-spectrum (FFT_B * conj(FFT_A))*/
static double       *xcorr     = NULL; /* Raw cross-correlation (time domain) */

static fftw_plan     pa = NULL;        /* FFTW plan: signal A -> FFT          */
static fftw_plan     pb = NULL;        /* FFTW plan: signal B -> FFT          */
static fftw_plan     px = NULL;        /* FFTW plan: cross-spectrum -> IFFT   */

/* -----------------------------------------------------------------------
 * Internal helpers for FFT buffer management
 * -----------------------------------------------------------------------*/

/** Free all dynamically allocated FFT buffers. */
static void _xcorr_free(void)
{
    if (xspec)    fftw_free(xspec);
    if (fftb)     fftw_free(fftb);
    if (ffta)     fftw_free(ffta);

    if (lags_loc) free(lags_loc);
    if (xcorr)    free(xcorr);
    if (sigb_loc) free(sigb_loc);
    if (siga_loc) free(siga_loc);

    xspec = NULL;   fftb = NULL;    ffta = NULL;
    lags_loc = NULL; xcorr = NULL;  sigb_loc = NULL; siga_loc = NULL;
}

/** Allocate all buffers needed for cross-correlation of length _N. */
static void _xcorr_malloc(void)
{
    siga_loc = calloc(_N, sizeof(double));
    sigb_loc = calloc(_N, sizeof(double));
    xcorr    = calloc(_N + 1, sizeof(double)); /* FFTW c2r needs N+1 */
    lags_loc = calloc(_N, sizeof(double));

    /* FFTW provides its own allocator that guarantees proper alignment
     * for SIMD instructions (SSE2, AVX, etc.).  Always use it for
     * fftw_complex arrays. */
    ffta  = fftw_alloc_complex(_N);
    fftb  = fftw_alloc_complex(_N);
    xspec = fftw_alloc_complex(_N);
}

/** Destroy existing FFTW plans and free all buffers. */
static void _xcorr_teardown(void)
{
    if (pa) fftw_destroy_plan(pa);
    if (pb) fftw_destroy_plan(pb);
    if (px) fftw_destroy_plan(px);
    pa = NULL; pb = NULL; px = NULL;

    _xcorr_free();
    fftw_cleanup();
}

/**
 * Allocate buffers and create FFTW plans for signals of length _N.
 *
 * FFTW plans describe *how* to compute an FFT of a given size.
 * Creating a plan is slow (FFTW tries different strategies), but
 * executing one is very fast.  That is why we cache them.
 *
 * We use real-to-complex (r2c) for forward transforms because our
 * input signals are real-valued, and complex-to-real (c2r) for the
 * inverse transform.
 */
static void _xcorr_setup(void)
{
    _xcorr_teardown();
    _xcorr_malloc();

    /* FFTW_ESTIMATE: pick a reasonable plan quickly (don't benchmark).
     * FFTW_DESTROY_INPUT: FFTW may overwrite the input array -- that's
     * fine because we always copy the data into our local buffers first. */
    pa = fftw_plan_dft_r2c_1d(_N, siga_loc, ffta,  FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    pb = fftw_plan_dft_r2c_1d(_N, sigb_loc, fftb,  FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    px = fftw_plan_dft_c2r_1d(_N, xspec,    xcorr, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
}

/**
 * Ensure the cached FFT setup matches the requested signal length.
 * If the length changed, rebuild everything.  Otherwise do nothing.
 */
static void _gcc_setup(unsigned N)
{
    if (_N != N) {
        _N = N;
        _xcorr_setup();
    }
}

/* -----------------------------------------------------------------------
 * Internal helpers for cross-correlation post-processing
 * -----------------------------------------------------------------------*/

/**
 * Shift an FFT output so that the zero-lag is in the middle.
 *
 * After an inverse FFT, the zero-lag (time=0) value is at index 0,
 * with positive lags at the start and negative lags at the end.
 * This function rearranges the array so that the zero-lag ends up
 * at index ceil(N/2), which is more intuitive for reading off delays.
 *
 * Before: [0, +1, +2, ..., +N/2, -N/2+1, ..., -2, -1]
 * After:  [-N/2+1, ..., -1, 0, +1, ..., +N/2]
 */
static void _fftshift(const double *in, double *out, unsigned N)
{
    /* The split point: how many elements are in the "second half" */
    unsigned half = (unsigned)floor((double)N / 2.0);

    /* Copy the second half of in[] (negative lags) to the start of out[] */
    memcpy(out, &in[half], sizeof(double) * (N - half));

    /* Copy the first half of in[] (positive lags) to the end of out[] */
    memcpy(&out[N - half], in, sizeof(double) * half);
}

/**
 * Find the index and value of the maximum element in an array.
 *
 * @param a     Input array.
 * @param N     Length of the array (must be >= 1).
 * @param max   Output: the maximum value found.
 * @param maxi  Output: the index of that maximum value.
 */
static void _max_index(const double *a, unsigned N,
                       double *max, unsigned *maxi)
{
    assert(a != NULL);
    assert(N >= 1);

    unsigned best_i = 0;
    double   best_v = a[0];

    for (unsigned i = 1; i < N; i++) {
        if (a[i] > best_v) {
            best_v = a[i];
            best_i = i;
        }
    }

    *max  = best_v;
    *maxi = best_i;
}

/* -----------------------------------------------------------------------
 * Public API: resource cleanup
 * -----------------------------------------------------------------------*/

void MD_shutdown(void)
{
    _xcorr_teardown();
}

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
    assert(a != NULL);
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
    assert(a != NULL);
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
    assert(a != NULL);
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
    assert(in != NULL);
    assert(out != NULL);
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
    assert(in != NULL);
    assert(out != NULL);
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
    assert(a != NULL);

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

/* -----------------------------------------------------------------------
 * Public API: GCC-PHAT delay estimation
 * -----------------------------------------------------------------------*/

/**
 * Estimate delays between a reference signal and several other signals.
 *
 * This is a convenience wrapper that calls MD_get_delay() for each
 * non-reference signal.  A Hanning window is applied to all signals
 * before computing the cross-correlation.
 *
 * @param sigs        Array of M pointers to signals.  sigs[0] is the reference.
 * @param M           Number of signals (must be >= 2).
 * @param N           Length of every signal.
 * @param margin      Search +/- margin samples around zero-lag.
 * @param weightfunc  SIMP or PHAT.
 * @param outdelays   Output: M-1 delay values (pre-allocated by caller).
 */
void MD_get_multiple_delays(const double **sigs, unsigned M, unsigned N,
                            unsigned margin, int weightfunc,
                            int *outdelays)
{
    /* Static buffers are reused across calls for efficiency.
     * They are only re-allocated when the signal length changes. */
    static unsigned last_N   = 0;
    static double  *hann_win = NULL;
    static double  *t_ref    = NULL;
    static double  *t_sig    = NULL;

    if (M < 2) return;

    /* Re-allocate if the signal length has changed */
    if (last_N != N) {
        free(hann_win);  /* free(NULL) is safe in C */
        free(t_ref);
        free(t_sig);

        hann_win = malloc(N * sizeof(double));
        assert(hann_win != NULL);
        MD_Gen_Hann_Win(hann_win, N);

        t_ref = malloc(N * sizeof(double));
        assert(t_ref != NULL);

        t_sig = malloc(N * sizeof(double));
        assert(t_sig != NULL);

        last_N = N;
    }

    /* Apply the Hanning window to the reference signal */
    for (unsigned i = 0; i < N; i++) {
        t_ref[i] = sigs[0][i] * hann_win[i];
    }

    /* Compute the delay for each non-reference signal */
    for (unsigned i = 0; i < M - 1; i++) {
        /* Window the comparison signal */
        for (unsigned j = 0; j < N; j++) {
            t_sig[j] = sigs[i + 1][j] * hann_win[j];
        }
        outdelays[i] = MD_get_delay(t_ref, t_sig, N, NULL, margin, weightfunc);
    }
}

/**
 * Estimate the delay (in samples) between two signals.
 *
 * The function computes the GCC between the two signals, then searches
 * within a +/- margin window around zero-lag for the peak.  The offset
 * of that peak from zero is the estimated delay.
 *
 * @param siga        First signal (reference).
 * @param sigb        Second signal.
 * @param N           Length of both signals.
 * @param ent         If non-NULL, receives the normalised entropy of the
 *                    lag values in the search window.  High entropy (~1.0)
 *                    means a flat, unreliable correlation; low entropy (~0.0)
 *                    means a sharp, trustworthy peak.
 * @param margin      How far (in samples) to search around zero-lag.
 * @param weightfunc  SIMP or PHAT.
 * @return            Delay in samples.
 */
int MD_get_delay(const double *siga, const double *sigb, unsigned N,
                 double *ent, unsigned margin, int weightfunc)
{
    _gcc_setup(N);

    /* Compute the full cross-correlation */
    MD_gcc(siga, sigb, N, lags_loc, weightfunc);

    /* The zero-lag position in the shifted output */
    unsigned center = (unsigned)ceil((double)N / 2.0);

    /* Clamp the margin so we don't read outside the array */
    unsigned m = margin;
    if (center < m) {
        m = center;
    }
    if (center + m >= N) {
        m = (N - 1) - center;
    }

    /* Search the window [center-m, center+m] for the peak */
    unsigned start = center - m;
    unsigned len   = 2 * m + 1;
    double   peak_val;
    unsigned peak_i;

    _max_index(lags_loc + start, len, &peak_val, &peak_i);

    /* Optionally compute the entropy of the search window */
    if (ent != NULL) {
        *ent = MD_entropy(lags_loc + start, len, true);
    }

    /* Convert the peak index to a signed delay relative to zero-lag */
    return (int)(peak_i) - (int)(m);
}

/**
 * Compute the Generalized Cross-Correlation between two signals.
 *
 * This is the core GCC-PHAT routine.  It fills lagvals[] with the
 * cross-correlation values, shifted so that the zero-lag value is
 * at index ceil(N/2).
 *
 * @param siga        First signal.
 * @param sigb        Second signal.
 * @param N           Length of both signals.
 * @param lagvals     Output array of N doubles (pre-allocated by caller).
 * @param weightfunc  SIMP or PHAT.
 */
void MD_gcc(const double *siga, const double *sigb, unsigned N,
            double *lagvals, int weightfunc)
{
    _gcc_setup(N);

    /* Copy input into local buffers (FFTW may overwrite them) */
    memcpy(siga_loc, siga, _N * sizeof(double));
    memcpy(sigb_loc, sigb, _N * sizeof(double));

    /* Forward FFT of both signals.
     *
     * For a real-valued input of length N, the FFT output has only
     * N/2 + 1 unique complex values (the rest are conjugate-symmetric).
     * FFTW's r2c transform exploits this and only writes N/2+1 values. */
    fftw_execute(pa);
    fftw_execute(pb);

    /* Compute the weighted cross-spectrum.
     *
     * The cross-spectrum is:  X[k] = FFT_B[k] * conj(FFT_A[k])
     *
     * For PHAT weighting, we divide by the magnitude of the cross-spectrum.
     * This keeps only the phase information, which produces a much sharper
     * peak in the time-domain correlation.
     *
     * We only loop over the N/2+1 non-redundant frequency bins. */
    unsigned num_bins = _N / 2 + 1;

    switch (weightfunc) {
    case PHAT:
        for (unsigned i = 0; i < num_bins; i++) {
            double complex cs = fftb[i] * conj(ffta[i]);
            /* Divide by magnitude; add DBL_MIN to prevent division by zero */
            xspec[i] = cs / (cabs(cs) + DBL_MIN);
        }
        break;

    default: /* SIMP: simple 1/N weighting */
        for (unsigned i = 0; i < num_bins; i++) {
            double complex cs = fftb[i] * conj(ffta[i]);
            xspec[i] = cs / (double)_N;
        }
        break;
    }

    /* Inverse FFT of the cross-spectrum -> time-domain cross-correlation */
    fftw_execute(px);

    /* Rearrange so the zero-lag is at the centre of the output array */
    _fftshift(xcorr, lagvals, _N);
}
