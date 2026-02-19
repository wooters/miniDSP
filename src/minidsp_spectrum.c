/**
 * @file minidsp_spectrum.c
 * @brief FFT-based spectrum analysis: magnitude spectrum, PSD, phase
 *        spectrum, STFT, mel filterbanks, MFCCs, FFT-based F0 estimation,
 *        and MD_shutdown().
 * @author Chuck Wooters <wooters@hey.com>
 * @copyright 2013 International Computer Science Institute
 */

#include "minidsp.h"
#include "minidsp_internal.h"

/* -----------------------------------------------------------------------
 * Static (file-scope) variables for spectrum analysis caching
 *
 * These are separate from the GCC cache (in minidsp_gcc.c) so that
 * spectrum analysis and delay estimation can be used independently
 * without invalidating each other's plans.
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

/* -----------------------------------------------------------------------
 * Static variables for mel filterbank cache
 *
 * MFCC extraction is often called frame-by-frame with fixed parameters.
 * Cache the triangular mel matrix to avoid rebuilding it on every frame.
 * -----------------------------------------------------------------------*/

static unsigned _mel_N = 0;
static unsigned _mel_num_mels = 0;
static double _mel_sample_rate = 0.0;
static double _mel_min_freq_hz = 0.0;
static double _mel_max_freq_hz = 0.0;
static double *_mel_fb = nullptr;   /* Row-major: [num_mels][N/2+1] */

/** F0-FFT peak must stand out from the average in-range spectrum. */
#define MD_F0_FFT_PROMINENCE_RATIO 2.5
/** Candidate must also be a meaningful fraction of the frame's global peak. */
#define MD_F0_FFT_GLOBAL_RATIO 0.2
/** Floor before logarithm in MFCC log-mel compression. */
#define MD_MFCC_LOG_FLOOR 1e-12

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

/** HTK mel scale conversion: Hz -> mel. */
static double md_hz_to_mel(double hz)
{
    return 2595.0 * log10(1.0 + hz / 700.0);
}

/** HTK mel scale conversion: mel -> Hz. */
static double md_mel_to_hz(double mel)
{
    return 700.0 * (pow(10.0, mel / 2595.0) - 1.0);
}

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

/** Free mel filterbank cache. */
static void _mel_teardown(void)
{
    free(_mel_fb);
    _mel_fb = nullptr;
    _mel_N = 0;
    _mel_num_mels = 0;
    _mel_sample_rate = 0.0;
    _mel_min_freq_hz = 0.0;
    _mel_max_freq_hz = 0.0;
}

/** Build the cached mel filterbank for the requested parameters. */
static void _mel_setup(unsigned N, double sample_rate, unsigned num_mels,
                       double min_freq_hz, double max_freq_hz)
{
    double nyquist = 0.5 * sample_rate;
    double f_lo = fmax(0.0, min_freq_hz);
    double f_hi = fmin(max_freq_hz, nyquist);

    if (_mel_fb != nullptr
        && _mel_N == N
        && _mel_num_mels == num_mels
        && _mel_sample_rate == sample_rate
        && _mel_min_freq_hz == f_lo
        && _mel_max_freq_hz == f_hi) {
        return;
    }

    _mel_teardown();

    _mel_N = N;
    _mel_num_mels = num_mels;
    _mel_sample_rate = sample_rate;
    _mel_min_freq_hz = f_lo;
    _mel_max_freq_hz = f_hi;

    unsigned num_bins = N / 2 + 1;
    size_t fb_len = (size_t)num_mels * num_bins;
    _mel_fb = calloc(fb_len, sizeof(double));
    assert(_mel_fb != nullptr);

    if (f_hi <= f_lo) return;

    unsigned num_points = num_mels + 2;
    double *mel_points = malloc(num_points * sizeof(double));
    double *hz_points = malloc(num_points * sizeof(double));
    assert(mel_points != nullptr);
    assert(hz_points != nullptr);

    double mel_min = md_hz_to_mel(f_lo);
    double mel_max = md_hz_to_mel(f_hi);
    if (mel_max <= mel_min) {
        free(hz_points);
        free(mel_points);
        return;
    }

    double mel_step = (mel_max - mel_min) / (double)(num_mels + 1);
    for (unsigned i = 0; i < num_points; i++) {
        mel_points[i] = mel_min + mel_step * (double)i;
        hz_points[i] = md_mel_to_hz(mel_points[i]);
    }

    for (unsigned m = 0; m < num_mels; m++) {
        double f_left = hz_points[m];
        double f_center = hz_points[m + 1];
        double f_right = hz_points[m + 2];
        if (f_center <= f_left || f_right <= f_center) {
            continue;
        }

        for (unsigned k = 0; k < num_bins; k++) {
            double f_bin = (double)k * sample_rate / (double)N;
            double w = 0.0;
            if (f_bin > f_left && f_bin < f_right) {
                if (f_bin <= f_center) {
                    w = (f_bin - f_left) / (f_center - f_left);
                } else {
                    w = (f_right - f_bin) / (f_right - f_center);
                }
            }
            if (w < 0.0) w = 0.0;
            if (w > 1.0) w = 1.0;
            _mel_fb[(size_t)m * num_bins + k] = w;
        }
    }

    free(hz_points);
    free(mel_points);
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
 * Public API: resource cleanup
 * -----------------------------------------------------------------------*/

void MD_shutdown(void)
{
    md_gcc_teardown();
    _mel_teardown();
    _stft_teardown();
    _spec_teardown();
    _spec_N = 0;
    fftw_cleanup();
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

double MD_f0_fft(const double *signal, unsigned N,
                 double sample_rate,
                 double min_freq_hz, double max_freq_hz)
{
    assert(signal != nullptr);
    assert(N >= 2);
    assert(sample_rate > 0.0);
    assert(min_freq_hz > 0.0);
    assert(max_freq_hz > min_freq_hz);

    unsigned k_min = (unsigned)ceil(min_freq_hz * (double)N / sample_rate);
    unsigned k_max = (unsigned)floor(max_freq_hz * (double)N / sample_rate);

    if (k_min < 1) k_min = 1;           /* skip DC */
    if (k_max > N / 2) k_max = N / 2;   /* one-sided Nyquist bound */
    if (k_min > k_max) return 0.0;

    if (MD_energy(signal, N) == 0.0) return 0.0;

    /* Reuse the shared STFT Hann cache to window this frame. */
    _stft_setup(N);
    for (unsigned n = 0; n < N; n++) {
        _spec_in[n] = signal[n] * _stft_win[n];
    }
    fftw_execute(_spec_plan);

    double global_max = 0.0;
    for (unsigned k = 1; k <= N / 2; k++) {
        double mag = cabs(_spec_out[k]);
        if (mag > global_max) global_max = mag;
    }

    double sum_mag = 0.0;
    double best_mag = -DBL_MAX;
    unsigned best_k = 0;
    unsigned num_bins = 0;

    for (unsigned k = k_min; k <= k_max; k++) {
        double mag = cabs(_spec_out[k]);
        sum_mag += mag;
        num_bins++;
        if (mag > best_mag) {
            best_mag = mag;
            best_k = k;
        }
    }

    if (best_k == 0 || best_mag <= 0.0 || num_bins == 0) return 0.0;

    double mean_mag = sum_mag / (double)num_bins;
    if (mean_mag <= 0.0) return 0.0;
    if (best_mag < MD_F0_FFT_PROMINENCE_RATIO * mean_mag) return 0.0;
    if (global_max <= 0.0) return 0.0;
    if (best_mag < MD_F0_FFT_GLOBAL_RATIO * global_max) return 0.0;

    double k_est = (double)best_k;
    if (best_k > 0 && best_k < N / 2) {
        double m_left  = cabs(_spec_out[best_k - 1]);
        double m_mid   = cabs(_spec_out[best_k]);
        double m_right = cabs(_spec_out[best_k + 1]);
        k_est += md_parabolic_offset(m_left, m_mid, m_right);
    }

    if (k_est < (double)k_min) k_est = (double)k_min;
    if (k_est > (double)k_max) k_est = (double)k_max;

    return k_est * sample_rate / (double)N;
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

void MD_mel_filterbank(unsigned N, double sample_rate,
                       unsigned num_mels,
                       double min_freq_hz, double max_freq_hz,
                       double *filterbank_out)
{
    assert(N >= 2);
    assert(sample_rate > 0.0);
    assert(num_mels > 0);
    assert(min_freq_hz < max_freq_hz);
    assert(filterbank_out != nullptr);

    _mel_setup(N, sample_rate, num_mels, min_freq_hz, max_freq_hz);

    unsigned num_bins = N / 2 + 1;
    size_t fb_len = (size_t)num_mels * num_bins;
    memcpy(filterbank_out, _mel_fb, fb_len * sizeof(double));
}

void MD_mel_energies(const double *signal, unsigned N,
                     double sample_rate, unsigned num_mels,
                     double min_freq_hz, double max_freq_hz,
                     double *mel_out)
{
    assert(signal != nullptr);
    assert(mel_out != nullptr);
    assert(N >= 2);
    assert(sample_rate > 0.0);
    assert(num_mels > 0);
    assert(min_freq_hz < max_freq_hz);

    _stft_setup(N);
    _mel_setup(N, sample_rate, num_mels, min_freq_hz, max_freq_hz);

    for (unsigned n = 0; n < N; n++) {
        _spec_in[n] = signal[n] * _stft_win[n];
    }
    fftw_execute(_spec_plan);

    unsigned num_bins = N / 2 + 1;
    for (unsigned m = 0; m < num_mels; m++) {
        mel_out[m] = 0.0;
    }

    for (unsigned k = 0; k < num_bins; k++) {
        double re = creal(_spec_out[k]);
        double im = cimag(_spec_out[k]);
        double psd = (re * re + im * im) / (double)N;

        for (unsigned m = 0; m < num_mels; m++) {
            mel_out[m] += _mel_fb[(size_t)m * num_bins + k] * psd;
        }
    }
}

void MD_mfcc(const double *signal, unsigned N,
             double sample_rate,
             unsigned num_mels, unsigned num_coeffs,
             double min_freq_hz, double max_freq_hz,
             double *mfcc_out)
{
    assert(signal != nullptr);
    assert(mfcc_out != nullptr);
    assert(N >= 2);
    assert(sample_rate > 0.0);
    assert(num_mels > 0);
    assert(num_coeffs >= 1 && num_coeffs <= num_mels);
    assert(min_freq_hz < max_freq_hz);

    double *mel = malloc(num_mels * sizeof(double));
    double *log_mel = malloc(num_mels * sizeof(double));
    assert(mel != nullptr);
    assert(log_mel != nullptr);

    MD_mel_energies(signal, N, sample_rate, num_mels,
                    min_freq_hz, max_freq_hz, mel);

    for (unsigned m = 0; m < num_mels; m++) {
        log_mel[m] = log(fmax(mel[m], MD_MFCC_LOG_FLOOR));
    }

    for (unsigned c = 0; c < num_coeffs; c++) {
        double norm = (c == 0)
                    ? sqrt(1.0 / (double)num_mels)
                    : sqrt(2.0 / (double)num_mels);
        double sum = 0.0;
        for (unsigned m = 0; m < num_mels; m++) {
            double angle = M_PI * (double)c * ((double)m + 0.5)
                         / (double)num_mels;
            sum += log_mel[m] * cos(angle);
        }
        mfcc_out[c] = norm * sum;
    }

    free(log_mel);
    free(mel);
}
