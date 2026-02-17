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
static double       *siga_loc  = nullptr; /* Local copy of input signal A        */
static double       *sigb_loc  = nullptr; /* Local copy of input signal B        */
static double       *lags_loc  = nullptr; /* Shifted cross-correlation output    */
static fftw_complex *ffta      = nullptr; /* FFT of signal A                     */
static fftw_complex *fftb      = nullptr; /* FFT of signal B                     */
static fftw_complex *xspec     = nullptr; /* Cross-spectrum (FFT_B * conj(FFT_A))*/
static double       *xcorr     = nullptr; /* Raw cross-correlation (time domain) */

static fftw_plan     pa = nullptr;        /* FFTW plan: signal A -> FFT          */
static fftw_plan     pb = nullptr;        /* FFTW plan: signal B -> FFT          */
static fftw_plan     px = nullptr;        /* FFTW plan: cross-spectrum -> IFFT   */

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

    xspec = nullptr;   fftb = nullptr;    ffta = nullptr;
    lags_loc = nullptr; xcorr = nullptr;  sigb_loc = nullptr; siga_loc = nullptr;
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
    pa = nullptr; pb = nullptr; px = nullptr;

    _xcorr_free();
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
 * Static (file-scope) variables for spectrum analysis caching
 *
 * These are separate from the GCC cache above so that spectrum analysis
 * and delay estimation can be used independently without invalidating
 * each other's plans.
 * -----------------------------------------------------------------------*/

static unsigned      _spec_N       = 0;       /* Cached signal length      */
static double       *_spec_in      = nullptr;  /* Local copy of input       */
static fftw_complex *_spec_out     = nullptr;  /* FFT output                */
static fftw_plan     _spec_plan    = nullptr;  /* FFTW r2c plan             */

/* -----------------------------------------------------------------------
 * Static variables for STFT Hanning window cache
 *
 * The STFT reuses the shared _spec_* r2c plan above.  Only the Hanning
 * window buffer is separate, since different callers may use different N.
 * -----------------------------------------------------------------------*/

static double  *_stft_win   = nullptr; /* Cached Hanning window             */
static unsigned _stft_win_N = 0;       /* Window length corresponding to above */

/** Free spectrum analysis buffers. */
static void _spec_free(void)
{
    if (_spec_out) fftw_free(_spec_out);
    if (_spec_in)  free(_spec_in);
    _spec_out = nullptr;
    _spec_in  = nullptr;
}

/** Tear down spectrum analysis plan and buffers. */
static void _spec_teardown(void)
{
    if (_spec_plan) fftw_destroy_plan(_spec_plan);
    _spec_plan = nullptr;
    _spec_free();
}

/**
 * Ensure the spectrum FFT cache matches the requested length.
 * Rebuilds the plan and buffers only when N changes.
 */
static void _spec_setup(unsigned N)
{
    if (_spec_N == N) return;

    _spec_teardown();
    _spec_N = N;

    _spec_in  = calloc(N, sizeof(double));
    _spec_out = fftw_alloc_complex(N / 2 + 1);

    _spec_plan = fftw_plan_dft_r2c_1d(
        (int)N, _spec_in, _spec_out,
        FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
}

/** Free the cached STFT Hanning window.  Called only from MD_shutdown(). */
static void _stft_teardown(void)
{
    free(_stft_win);
    _stft_win   = nullptr;
    _stft_win_N = 0;
}

/**
 * Ensure the shared r2c plan is ready for length N, and that the STFT
 * Hanning window buffer matches N.  Rebuilds the window only when N changes.
 */
static void _stft_setup(unsigned N)
{
    _spec_setup(N);                 /* reuse shared r2c plan */
    if (_stft_win_N == N) return;   /* fast path: window already correct */

    free(_stft_win);
    _stft_win = malloc(N * sizeof(double));
    assert(_stft_win != nullptr);
    MD_Gen_Hann_Win(_stft_win, N);
    _stft_win_N = N;
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
    assert(a != nullptr);
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
    _stft_teardown();
    _spec_teardown();
    _spec_N = 0;
    fftw_cleanup();
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

/* -----------------------------------------------------------------------
 * Public API: Signal generators
 * -----------------------------------------------------------------------*/

void MD_sine_wave(double *output, unsigned N, double amplitude,
                  double freq, double sample_rate)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);
    double phase_step = 2.0 * M_PI * freq / sample_rate;
    for (unsigned i = 0; i < N; i++)
        output[i] = amplitude * sin(phase_step * (double)i);
}

/* -----------------------------------------------------------------------
 * Public API: FFT / Spectrum Analysis
 * -----------------------------------------------------------------------*/

/**
 * Compute the magnitude spectrum of a real-valued signal.
 *
 * This performs a real-to-complex FFT using FFTW, then computes the
 * absolute value (magnitude) of each complex frequency bin.
 *
 * For a real signal of length N, the FFT is conjugate-symmetric, so
 * only the first N/2 + 1 bins are unique:
 *
 *   - Bin 0:     DC component (zero frequency)
 *   - Bin k:     frequency = k * sample_rate / N
 *   - Bin N/2:   Nyquist frequency (sample_rate / 2)
 *
 * The magnitude is computed as:
 *   |X(k)| = sqrt( Re(X(k))^2 + Im(X(k))^2 )
 *
 * The FFT plan is cached and reused across calls of the same length,
 * following the same pattern as the GCC functions.
 *
 * @param signal   Input signal of length N.
 * @param N        Number of samples (must be >= 2).
 * @param mag_out  Output: magnitudes for bins 0..N/2.
 *                 Must be pre-allocated to N/2 + 1 doubles.
 */
void MD_magnitude_spectrum(const double *signal, unsigned N, double *mag_out)
{
    assert(signal != nullptr);
    assert(mag_out != nullptr);
    assert(N >= 2);

    _spec_setup(N);

    /* Copy input into the local buffer (FFTW may overwrite it) */
    memcpy(_spec_in, signal, N * sizeof(double));

    /* Execute the forward FFT (real -> complex) */
    fftw_execute(_spec_plan);

    /* Compute magnitude |X(k)| = sqrt(re^2 + im^2) for each bin */
    unsigned num_bins = N / 2 + 1;
    for (unsigned k = 0; k < num_bins; k++) {
        mag_out[k] = cabs(_spec_out[k]);
    }
}

/**
 * Compute the power spectral density (PSD) of a real-valued signal.
 *
 * The PSD is the "periodogram" estimator: PSD[k] = |X(k)|^2 / N.
 * It reuses the same FFT cache as MD_magnitude_spectrum() -- both
 * perform the same real-to-complex FFT, only the post-processing differs.
 *
 * We compute |X(k)|^2 = re^2 + im^2 directly from the real and imaginary
 * parts, rather than calling cabs() (which computes sqrt(re^2 + im^2))
 * and then squaring.  This avoids a redundant sqrt and is both faster
 * and more numerically precise.
 *
 * @param signal   Input signal of length N.
 * @param N        Number of samples (must be >= 2).
 * @param psd_out  Output: PSD for bins 0..N/2.
 *                 Must be pre-allocated to N/2 + 1 doubles.
 */
void MD_power_spectral_density(const double *signal, unsigned N, double *psd_out)
{
    assert(signal != nullptr);
    assert(psd_out != nullptr);
    assert(N >= 2);

    _spec_setup(N);

    /* Copy input into the local buffer (FFTW may overwrite it) */
    memcpy(_spec_in, signal, N * sizeof(double));

    /* Execute the forward FFT (real -> complex) */
    fftw_execute(_spec_plan);

    /* Compute PSD[k] = |X(k)|^2 / N = (re^2 + im^2) / N.
     * We use creal/cimag instead of cabs to avoid a redundant sqrt --
     * cabs computes sqrt(re^2 + im^2), and we'd just square it again. */
    unsigned num_bins = N / 2 + 1;
    for (unsigned k = 0; k < num_bins; k++) {
        double re = creal(_spec_out[k]);
        double im = cimag(_spec_out[k]);
        psd_out[k] = (re * re + im * im) / (double)N;
    }
}

/**
 * Compute the one-sided phase spectrum of a real-valued signal.
 *
 * The phase is the argument (angle) of each complex DFT coefficient:
 *   phi(k) = atan2(Im(X(k)), Re(X(k)))
 *
 * This reuses the same FFT cache as MD_magnitude_spectrum() and
 * MD_power_spectral_density() -- no additional plan allocation.
 *
 * Phase is scale-invariant (multiplying a signal by a positive constant
 * does not change its phase), so no normalisation by N is needed.
 *
 * @param signal    Input signal of length N.
 * @param N         Number of samples (must be >= 2).
 * @param phase_out Output: phase in radians for bins 0..N/2.
 *                  Must be pre-allocated to N/2 + 1 doubles.
 */
void MD_phase_spectrum(const double *signal, unsigned N, double *phase_out)
{
    assert(signal    != nullptr);
    assert(phase_out != nullptr);
    assert(N >= 2);

    _spec_setup(N);

    /* Copy input into the local buffer (FFTW may overwrite it) */
    memcpy(_spec_in, signal, N * sizeof(double));

    /* Execute the forward FFT (real -> complex) */
    fftw_execute(_spec_plan);

    /* Compute phi(k) = atan2(Im(X(k)), Re(X(k))) for each bin.
     * carg() from <complex.h> does exactly this and returns [-pi, pi]. */
    unsigned num_bins = N / 2 + 1;
    for (unsigned k = 0; k < num_bins; k++) {
        phase_out[k] = carg(_spec_out[k]);
    }
}

/**
 * Compute the number of complete STFT frames for the given parameters.
 *
 * Returns (signal_len - N) / hop + 1 when signal_len >= N, else 0.
 * Integer division truncates, so only complete frames are counted.
 */
unsigned MD_stft_num_frames(unsigned signal_len, unsigned N, unsigned hop)
{
    if (signal_len < N) return 0;
    return (signal_len - N) / hop + 1;
}

/**
 * Compute the Short-Time Fourier Transform (STFT) magnitude matrix.
 *
 * Slides a Hanning-windowed r2c FFT over the signal in steps of @p hop
 * samples.  For each frame, the window is applied by multiplying directly
 * into the shared _spec_in buffer (no separate memcpy pass), then FFTW
 * executes the plan and magnitudes |X(k)| are written to the output row.
 *
 * The shared _spec_* plan is reused so that interleaving calls to
 * MD_stft(), MD_magnitude_spectrum(), and MD_power_spectral_density()
 * with the same N incurs no extra plan-rebuild overhead.
 *
 * @param signal      Input signal.
 * @param signal_len  Number of samples (0 frames if < N, see header).
 * @param N           FFT window size (>= 2).
 * @param hop         Hop size (>= 1).
 * @param mag_out     Row-major output: mag_out[f*(N/2+1) + k] = |X_f(k)|.
 *                    Must be pre-allocated by caller.
 */
void MD_stft(const double *signal, unsigned signal_len,
             unsigned N, unsigned hop, double *mag_out)
{
    assert(signal  != nullptr);
    assert(mag_out != nullptr);
    assert(N   >= 2);
    assert(hop >= 1);

    unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
    if (num_frames == 0) return;   /* signal too short for even one frame */

    _stft_setup(N);
    unsigned num_bins = N / 2 + 1;

    for (unsigned f = 0; f < num_frames; f++) {
        const double *src = signal + (size_t)f * hop;  /* (size_t) avoids overflow */

        /* Apply Hanning window directly into the shared input buffer.
         * Touching every element anyway, so no separate memcpy is needed. */
        for (unsigned n = 0; n < N; n++) {
            _spec_in[n] = src[n] * _stft_win[n];
        }

        fftw_execute(_spec_plan);

        /* Write magnitude row: |X_f(k)| = sqrt(re^2 + im^2).
         * cabs() is correct here -- we need the actual sqrt magnitude,
         * not |z|^2, so the creal/cimag shortcut from PSD does not apply. */
        double *row = mag_out + (size_t)f * num_bins;
        for (unsigned k = 0; k < num_bins; k++) {
            row[k] = cabs(_spec_out[k]);
        }
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
    static double  *hann_win = nullptr;
    static double  *t_ref    = nullptr;
    static double  *t_sig    = nullptr;

    if (M < 2) return;

    /* Re-allocate if the signal length has changed */
    if (last_N != N) {
        free(hann_win);  /* free(nullptr) is safe in C */
        free(t_ref);
        free(t_sig);

        hann_win = malloc(N * sizeof(double));
        assert(hann_win != nullptr);
        MD_Gen_Hann_Win(hann_win, N);

        t_ref = malloc(N * sizeof(double));
        assert(t_ref != nullptr);

        t_sig = malloc(N * sizeof(double));
        assert(t_sig != nullptr);

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
        outdelays[i] = MD_get_delay(t_ref, t_sig, N, nullptr, margin, weightfunc);
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
 * @param ent         If non-null, receives the normalised entropy of the
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
    if (ent != nullptr) {
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
