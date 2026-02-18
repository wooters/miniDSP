/**
 * @file minidsp.h
 * @brief A mini library of DSP (Digital Signal Processing) routines.
 *
 * This header declares functions for:
 *   - Basic signal measurements (energy, power, entropy)
 *   - Signal analysis (RMS, zero-crossing rate, autocorrelation, peak detection, mixing)
 *   - Signal scaling and gain adjustment
 *   - Window generation (Hanning, Hamming, Blackman, rectangular)
 *   - Signal generators (sine, white noise, impulse, chirp, square, sawtooth)
 *   - FIR filtering and convolution (time-domain and FFT overlap-add)
 *   - FFT-based magnitude spectrum, power spectral density, and STFT
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
#include <math.h>

/* M_PI is not guaranteed by the C standard.  Define it if the system didn't. */
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
 * Signal analysis
 * -----------------------------------------------------------------------*/

/**
 * Compute the root mean square (RMS) of a signal.
 *
 * RMS is the standard measure of signal "loudness":
 *
 * \f[
 * \mathrm{RMS} = \sqrt{\frac{1}{N}\sum_{n=0}^{N-1} x[n]^2}
 * \f]
 *
 * Equivalently, `sqrt(MD_power(a, N))`.  A unit-amplitude sine wave
 * has RMS = 1/sqrt(2) ~ 0.707.  A DC signal of value c has RMS = |c|.
 *
 * @param a  Input signal of length N.
 * @param N  Number of samples.  Must be > 0.
 * @return   RMS value (always >= 0).
 *
 * @code
 * double sig[1024];
 * MD_sine_wave(sig, 1024, 1.0, 440.0, 44100.0);
 * double rms = MD_rms(sig, 1024);  // ~ 0.707
 * @endcode
 */
double MD_rms(const double *a, unsigned N);

/**
 * Compute the zero-crossing rate of a signal.
 *
 * The zero-crossing rate counts how often the signal changes sign,
 * normalised by the number of adjacent pairs:
 *
 * \f[
 * \mathrm{ZCR} = \frac{1}{N-1}\sum_{n=1}^{N-1}
 *     \mathbf{1}\!\bigl[\mathrm{sgn}(x[n]) \ne \mathrm{sgn}(x[n-1])\bigr]
 * \f]
 *
 * Returns a value in [0.0, 1.0].  High ZCR indicates noise or
 * high-frequency content; low ZCR indicates tonal or low-frequency
 * content.  Zero is treated as non-negative (standard convention).
 *
 * @param a  Input signal of length N.
 * @param N  Number of samples.  Must be > 1.
 * @return   Zero-crossing rate in [0.0, 1.0].
 *
 * @code
 * double sig[4096];
 * MD_sine_wave(sig, 4096, 1.0, 1000.0, 16000.0);
 * double zcr = MD_zero_crossing_rate(sig, 4096);
 * // zcr ~ 2 * 1000 / 16000 = 0.125
 * @endcode
 */
double MD_zero_crossing_rate(const double *a, unsigned N);

/**
 * Compute the normalised autocorrelation of a signal.
 *
 * The autocorrelation measures the similarity between a signal and a
 * delayed copy of itself.  This function computes the normalised
 * (biased) autocorrelation for lags 0 through max_lag-1:
 *
 * \f[
 * R[\tau] = \frac{1}{R[0]} \sum_{n=0}^{N-1-\tau} x[n]\,x[n+\tau]
 * \f]
 *
 * where \f$R[0] = \sum x[n]^2\f$ is the signal energy.
 * The output satisfies out[0] = 1.0 and |out[tau]| <= 1.0.
 * A silent signal (all zeros) produces all-zero output.
 *
 * @param a        Input signal of length N.
 * @param N        Number of samples.  Must be > 0.
 * @param out      Output array of length max_lag (caller-allocated).
 * @param max_lag  Number of lag values to compute.  Must be > 0 and < N.
 *
 * @code
 * double sig[1024];
 * MD_sine_wave(sig, 1024, 1.0, 100.0, 1000.0);
 * double acf[50];
 * MD_autocorrelation(sig, 1024, acf, 50);
 * // acf[0] = 1.0, acf[10] ~ 1.0 (lag = one period of 100 Hz at 1 kHz)
 * @endcode
 */
void MD_autocorrelation(const double *a, unsigned N,
                        double *out, unsigned max_lag);

/**
 * Detect peaks (local maxima) in a signal.
 *
 * A sample a[i] is a peak if it is strictly greater than both
 * immediate neighbours and above the given threshold.  The
 * min_distance parameter suppresses nearby secondary peaks:
 * after accepting a peak at index i, any candidate within
 * min_distance indices is skipped.
 *
 * Peaks are found left-to-right (deterministic tie-breaking).
 * Endpoint samples (i=0, i=N-1) are never peaks because they
 * lack two neighbours.
 *
 * @param a             Input signal of length N.
 * @param N             Number of samples.
 * @param threshold     Minimum value for a peak.
 * @param min_distance  Minimum index gap between accepted peaks (>= 1).
 * @param peaks_out     Output array of peak indices (caller-allocated,
 *                      worst case N elements).
 * @param num_peaks_out Receives the number of peaks found.
 *
 * @code
 * double sig[] = {0, 1, 3, 1, 0, 2, 5, 2, 0};
 * unsigned peaks[9], n;
 * MD_peak_detect(sig, 9, 0.0, 1, peaks, &n);
 * // n = 2, peaks = {2, 6}  (values 3 and 5)
 * @endcode
 */
void MD_peak_detect(const double *a, unsigned N, double threshold,
                    unsigned min_distance, unsigned *peaks_out,
                    unsigned *num_peaks_out);

/**
 * Mix (weighted sum) two signals.
 *
 * Computes the element-wise weighted sum:
 *
 * \f[
 * \mathrm{out}[n] = w_a \cdot a[n] + w_b \cdot b[n]
 * \f]
 *
 * In-place safe: the output buffer may alias either input.
 *
 * @param a    First input signal of length N.
 * @param b    Second input signal of length N.
 * @param out  Output signal of length N (caller-allocated).
 * @param N    Number of samples.
 * @param w_a  Weight for signal a.
 * @param w_b  Weight for signal b.
 *
 * @code
 * double sine[1024], noise[1024], mix[1024];
 * MD_sine_wave(sine, 1024, 1.0, 440.0, 44100.0);
 * MD_white_noise(noise, 1024, 0.1, 42);
 * MD_mix(sine, noise, mix, 1024, 0.8, 0.2);
 * @endcode
 */
void MD_mix(const double *a, const double *b, double *out,
            unsigned N, double w_a, double w_b);

/* -----------------------------------------------------------------------
 * FIR filters / convolution
 * -----------------------------------------------------------------------*/

/**
 * Compute the output length of a full linear convolution.
 *
 * For input length N and kernel length M, full convolution length is N+M-1.
 */
unsigned MD_convolution_num_samples(unsigned signal_len, unsigned kernel_len);

/**
 * Time-domain full linear convolution (direct sum-of-products).
 *
 * Computes:
 *   out[n] = sum_{k=0}^{kernel_len-1} signal[n-k] * kernel[k]
 * with out-of-range signal samples treated as zero.
 *
 * @param signal      Input signal of length signal_len.
 * @param signal_len  Number of input samples (must be > 0).
 * @param kernel      FIR kernel of length kernel_len.
 * @param kernel_len  Number of FIR taps (must be > 0).
 * @param out         Output buffer of length signal_len + kernel_len - 1.
 */
void MD_convolution_time(const double *signal, unsigned signal_len,
                         const double *kernel, unsigned kernel_len,
                         double *out);

/**
 * Causal moving-average FIR filter with zero-padded startup.
 *
 * Computes:
 *   out[n] = (1/window_len) * sum_{k=0}^{window_len-1} signal[n-k]
 * where out-of-range samples (n-k < 0) are treated as zero.
 *
 * @param signal      Input signal of length signal_len.
 * @param signal_len  Number of input samples (must be > 0).
 * @param window_len  Moving-average window length (must be > 0).
 * @param out         Output buffer of length signal_len.
 */
void MD_moving_average(const double *signal, unsigned signal_len,
                       unsigned window_len, double *out);

/**
 * Apply a causal FIR filter with arbitrary coefficients.
 *
 * Computes:
 *   out[n] = sum_{k=0}^{num_taps-1} coeffs[k] * signal[n-k]
 * with out-of-range signal samples treated as zero.
 *
 * @param signal      Input signal of length signal_len.
 * @param signal_len  Number of input samples (must be > 0).
 * @param coeffs      FIR coefficients of length num_taps.
 * @param num_taps    Number of FIR taps (must be > 0).
 * @param out         Output buffer of length signal_len.
 */
void MD_fir_filter(const double *signal, unsigned signal_len,
                   const double *coeffs, unsigned num_taps,
                   double *out);

/**
 * Full linear convolution using FFT overlap-add (offline).
 *
 * Produces the same output as MD_convolution_time() but is faster for
 * longer kernels by processing blocks in the frequency domain.
 *
 * @param signal      Input signal of length signal_len.
 * @param signal_len  Number of input samples (must be > 0).
 * @param kernel      FIR kernel of length kernel_len.
 * @param kernel_len  Number of FIR taps (must be > 0).
 * @param out         Output buffer of length signal_len + kernel_len - 1.
 */
void MD_convolution_fft_ola(const double *signal, unsigned signal_len,
                            const double *kernel, unsigned kernel_len,
                            double *out);

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
 * FFT / Spectrum Analysis
 * -----------------------------------------------------------------------*/

/**
 * Compute the magnitude spectrum of a real-valued signal.
 *
 * Given a signal of length N, this function computes the FFT and returns
 * the magnitude |X(k)| for each frequency bin.  Because the input is
 * real-valued, the FFT output is conjugate-symmetric, so only the first
 * N/2 + 1 bins are unique (from DC through Nyquist).
 *
 * To convert bin index k to a frequency in Hz:
 *   freq_k = k * sample_rate / N
 *
 * The output is **not** normalised by N -- divide each value by N to get
 * the "standard" DFT magnitude, or by N/2 (except DC and Nyquist) for
 * single-sided amplitude.
 *
 * @param signal   Input signal of length N.
 * @param N        Number of samples in the signal.
 * @param mag_out  Output array, must be pre-allocated to at least N/2 + 1
 *                 doubles.  On return, mag_out[k] = |X(k)|.
 *
 * @note The caller must allocate mag_out.  The required size is
 *       (N / 2 + 1) * sizeof(double).
 *
 * Example:
 * @code
 *   double signal[1024];
 *   // ... fill signal with audio samples ...
 *
 *   unsigned num_bins = 1024 / 2 + 1;  // = 513
 *   double *mag = malloc(num_bins * sizeof(double));
 *   MD_magnitude_spectrum(signal, 1024, mag);
 *
 *   // mag[0]   = DC component magnitude
 *   // mag[k]   = magnitude at frequency k * sample_rate / 1024
 *   // mag[512] = Nyquist frequency magnitude
 *
 *   free(mag);
 *   MD_shutdown();  // free cached FFT plans when done
 * @endcode
 */
void MD_magnitude_spectrum(const double *signal, unsigned N, double *mag_out);

/**
 * Compute the power spectral density (PSD) of a real-valued signal.
 *
 * The PSD describes how a signal's power is distributed across frequencies.
 * While the magnitude spectrum tells you the *amplitude* at each frequency,
 * the PSD tells you the *power* -- useful for noise analysis, SNR estimation,
 * and comparing signals of different lengths.
 *
 * This function computes the "periodogram" estimator:
 *
 *   PSD[k] = |X(k)|^2 / N  =  (Re(X(k))^2 + Im(X(k))^2) / N
 *
 * where X(k) is the DFT of the input signal (unnormalised, as computed by FFTW).
 *
 * **Relationship to the magnitude spectrum:**
 *   PSD[k] = |X(k)|^2 / N  =  (magnitude[k])^2 / N
 *
 * **Parseval's theorem (energy conservation):**
 *   The one-sided PSD sums to the total signal energy:
 *   PSD[0] + 2 * sum(PSD[1..N/2-1]) + PSD[N/2] = sum(x[n]^2)
 *
 * @param signal   Input signal of length N.
 * @param N        Number of samples in the signal (must be >= 2).
 * @param psd_out  Output array, must be pre-allocated to at least N/2 + 1
 *                 doubles.  On return, psd_out[k] = |X(k)|^2 / N.
 *
 * @note The caller must allocate psd_out.  The required size is
 *       (N / 2 + 1) * sizeof(double).
 *
 * Example:
 * @code
 *   double signal[1024];
 *   // ... fill signal with audio samples ...
 *
 *   unsigned num_bins = 1024 / 2 + 1;  // = 513
 *   double *psd = malloc(num_bins * sizeof(double));
 *   MD_power_spectral_density(signal, 1024, psd);
 *
 *   // psd[0]   = DC power
 *   // psd[k]   = power at frequency k * sample_rate / 1024
 *   // psd[512] = Nyquist frequency power
 *
 *   free(psd);
 *   MD_shutdown();  // free cached FFT plans when done
 * @endcode
 */
void MD_power_spectral_density(const double *signal, unsigned N, double *psd_out);

/**
 * @brief Compute the one-sided phase spectrum of a real signal.
 *
 * Returns the instantaneous phase angle \f$\phi(k) = \arg X(k)\f$ for each
 * DFT bin using an unnormalised real-to-complex FFT (FFTW r2c).  The phase
 * is expressed in **radians** in the range \f$[-\pi,\,\pi]\f$.
 *
 * For a real signal of length \f$N\f$, only the first \f$N/2+1\f$ bins carry
 * unique information (bins \f$N/2+1\ldots N-1\f$ are conjugate-symmetric
 * mirrors).  Accordingly, @p phase_out must be pre-allocated to hold at least
 * \f$N/2+1\f$ doubles.
 *
 * **Interpretation:**
 * - A pure cosine at an integer bin \f$k_0\f$ (exact period in \f$N\f$
 *   samples) produces \f$\phi(k_0) = 0\f$.
 * - A pure sine at the same bin produces \f$\phi(k_0) = -\pi/2\f$.
 * - A time-delayed signal exhibits **linear phase**: \f$\phi(k) = -2\pi k d/N\f$,
 *   where \f$d\f$ is the delay in samples.
 *
 * @note Phase values at bins where the magnitude is near zero are numerically
 *       unreliable.  Always examine MD_magnitude_spectrum() alongside the
 *       phase spectrum to identify significant bins.
 *
 * @param[in]  signal    Input signal of length @p N.
 * @param[in]  N         Signal length (must be >= 2).
 * @param[out] phase_out Pre-allocated array of at least @p N/2+1 doubles
 *                       that receives the phase in radians.
 *
 * **Example**
 * @code
 * unsigned N = 1024;
 * double *sig = malloc(N * sizeof(double));
 * // ... fill sig ...
 * unsigned num_bins = N / 2 + 1;
 * double *phase = malloc(num_bins * sizeof(double));
 * MD_phase_spectrum(sig, N, phase);
 * // phase[k] in [-M_PI, M_PI]
 * free(sig);
 * free(phase);
 * MD_shutdown();
 * @endcode
 *
 * @see MD_magnitude_spectrum(), MD_power_spectral_density(), MD_shutdown()
 */
void MD_phase_spectrum(const double *signal, unsigned N, double *phase_out);

/**
 * Compute the number of STFT frames for the given signal length and parameters.
 *
 * The formula is:
 *   num_frames = (signal_len - N) / hop + 1  when signal_len >= N
 *   num_frames = 0                            when signal_len < N
 *
 * Use this function to size the output buffer before calling MD_stft().
 *
 * @param signal_len  Total number of samples in the signal.
 * @param N           FFT window size (samples per frame).
 * @param hop         Hop size (samples between successive frame starts).
 * @return            Number of complete frames that fit in the signal.
 *
 * Example:
 * @code
 *   // 1 s of audio at 16 kHz, 32 ms window, 8 ms hop
 *   unsigned signal_len = 16000;
 *   unsigned N   = 512;   // 32 ms at 16 kHz
 *   unsigned hop = 128;   // 8 ms (75% overlap)
 *   unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
 *   // num_frames = (16000 - 512) / 128 + 1 = 121
 * @endcode
 */
unsigned MD_stft_num_frames(unsigned signal_len, unsigned N, unsigned hop);

/**
 * Compute the Short-Time Fourier Transform (STFT) of a real-valued signal.
 *
 * The STFT slides a Hanning-windowed FFT over the signal in steps of
 * @p hop samples, producing a time-frequency magnitude matrix.
 *
 * For each frame @p f starting at sample @p f * hop, the function:
 *   1. Multiplies the frame by a Hanning window.
 *   2. Computes the real-to-complex FFT using FFTW.
 *   3. Stores the magnitudes |X(k)| for bins 0..N/2 in the output.
 *
 * The STFT formula for frame @p f and bin @p k is:
 *
 *   X_f(k) = SUM_{n=0}^{N-1}  w[n] * x[f*hop + n] * e^{-j2pi*k*n/N}
 *
 *   mag_out[f * (N/2+1) + k] = |X_f(k)|
 *
 * where @p w[n] is the Hanning window.
 *
 * The output is **not** normalised by N -- divide each value by N to get
 * the "standard" DFT magnitude, consistent with MD_magnitude_spectrum().
 *
 * The STFT reuses the same cached FFT plan as MD_magnitude_spectrum() and
 * MD_power_spectral_density().  Only the Hanning window buffer is separate.
 *
 * @param signal      Input signal.
 * @param signal_len  Total number of samples in the signal.
 * @param N           FFT window size (must be >= 2).
 * @param hop         Hop size in samples (must be >= 1).
 * @param mag_out     Output array (row-major).  Must be pre-allocated to at
 *                    least MD_stft_num_frames(signal_len, N, hop) * (N/2+1)
 *                    doubles.  On return, mag_out[f*(N/2+1) + k] = |X_f(k)|.
 *
 * @note  If signal_len < N, the function returns immediately without writing
 *        any output (zero frames fit).  Use MD_stft_num_frames() to check
 *        in advance.
 *
 * @note  The caller must allocate mag_out.  A typical pattern:
 * @code
 *   unsigned N   = 512;
 *   unsigned hop = 128;
 *   unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);
 *   unsigned num_bins   = N / 2 + 1;
 *
 *   double *mag_out = malloc(num_frames * num_bins * sizeof(double));
 *   MD_stft(signal, signal_len, N, hop, mag_out);
 *
 *   // mag_out[f * num_bins + k] = |X_f(k)|
 *   // Convert bin k to Hz: freq_hz = k * sample_rate / N
 *   // Convert frame f to seconds: time_s = (double)(f * hop) / sample_rate
 *
 *   free(mag_out);
 *   MD_shutdown();  // free cached FFT plans when done
 * @endcode
 */
void MD_stft(const double *signal, unsigned signal_len,
             unsigned N, unsigned hop,
             double *mag_out);

/* -----------------------------------------------------------------------
 * Window generation
 * -----------------------------------------------------------------------*/

/**
 * Generate a Hanning (Hann) window of length n.
 * Tapers to zero at both ends and is a common default for FFT analysis.
 */
void MD_Gen_Hann_Win(double *out, unsigned n);

/**
 * Generate a Hamming window of length n.
 * Similar to Hanning, but with a slightly lower first sidelobe.
 */
void MD_Gen_Hamming_Win(double *out, unsigned n);

/**
 * Generate a Blackman window of length n.
 * Much lower sidelobes than Hanning/Hamming, with a wider main lobe.
 */
void MD_Gen_Blackman_Win(double *out, unsigned n);

/**
 * Generate a rectangular window of length n (all ones).
 * Useful as a baseline reference (equivalent to no tapering).
 */
void MD_Gen_Rect_Win(double *out, unsigned n);

/* -----------------------------------------------------------------------
 * Signal generators
 * -----------------------------------------------------------------------*/

/**
 * Generate a sine wave.
 *
 * Fills `output[i] = amplitude * sin(2π * freq * i / sample_rate)`
 * for i in [0, N).
 *
 * This is the simplest test signal in DSP — a pure tone at a single
 * frequency.  Use it to verify filter responses, check FFT bin
 * alignment, or provide a clean input for any processing chain.
 *
 * @param output      Output buffer of length N (caller-allocated).
 * @param N           Number of samples to generate.  Must be > 0.
 * @param amplitude   Peak amplitude of the sine wave (e.g. 1.0).
 * @param freq        Frequency in Hz.
 * @param sample_rate Sampling rate in Hz.  Must be > 0.
 *
 * @code
 * double sig[1024];
 * MD_sine_wave(sig, 1024, 1.0, 440.0, 44100.0);  // 440 Hz A4 tone
 * @endcode
 */
void MD_sine_wave(double *output, unsigned N, double amplitude,
                  double freq, double sample_rate);

/**
 * Generate Gaussian white noise.
 *
 * Fills `output` with normally distributed random samples (mean 0,
 * standard deviation `amplitude`) using the Box-Muller transform.
 * White noise has equal energy at all frequencies -- its power
 * spectral density is approximately flat.
 *
 * Use white noise to test filters, measure impulse responses, or
 * as an additive noise source for SNR experiments.
 *
 * @param output    Output buffer of length N (caller-allocated).
 * @param N         Number of samples to generate.  Must be > 0.
 * @param amplitude Standard deviation of the noise (e.g. 1.0).
 * @param seed      Seed for the random number generator.  Using the
 *                  same seed produces the same output sequence, which
 *                  is useful for reproducible tests.
 *
 * @code
 * double noise[4096];
 * MD_white_noise(noise, 4096, 1.0, 42);  // reproducible Gaussian noise
 * @endcode
 */
void MD_white_noise(double *output, unsigned N, double amplitude,
                    unsigned seed);

/**
 * Generate a discrete impulse (Kronecker delta).
 *
 * Fills the output buffer with zeros except at @p position, where the
 * value is set to @p amplitude.  A unit impulse (amplitude 1.0 at
 * position 0) is the identity element of convolution and has a
 * perfectly flat magnitude spectrum.
 *
 * Common uses:
 *   - Measure a system's impulse response by feeding it through a filter.
 *   - Verify that MD_magnitude_spectrum() returns a flat spectrum.
 *   - Create delayed spikes for testing convolution and delay estimation.
 *
 * @param output    Output buffer of length N (caller-allocated).
 * @param N         Number of samples to generate.  Must be > 0.
 * @param amplitude Value of the impulse spike (e.g. 1.0 for unit impulse).
 * @param position  Sample index of the spike (0-based).  Must be < N.
 *
 * @code
 * double sig[1024];
 * MD_impulse(sig, 1024, 1.0, 0);  // unit impulse at sample 0
 * @endcode
 */
void MD_impulse(double *output, unsigned N, double amplitude, unsigned position);

/**
 * Generate a linear chirp (swept sine with linearly increasing frequency).
 *
 * The instantaneous frequency sweeps linearly from @p f_start to @p f_end
 * over N samples:
 *
 *   f(t) = f_start + (f_end - f_start) * t / T
 *
 * where T = (N-1) / sample_rate is the sweep duration.  The output is:
 *
 *   output[i] = amplitude * sin(2*pi * (f_start*t + 0.5*chirp_rate*t^2))
 *
 * A linear chirp is the standard test signal for spectrograms -- its
 * instantaneous frequency traces a straight diagonal line in the
 * time-frequency plane.
 *
 * @param output      Output buffer of length N (caller-allocated).
 * @param N           Number of samples to generate.  Must be > 0.
 * @param amplitude   Peak amplitude of the chirp (e.g. 1.0).
 * @param f_start     Starting frequency in Hz.
 * @param f_end       Ending frequency in Hz.
 * @param sample_rate Sampling rate in Hz.  Must be > 0.
 *
 * @code
 * double sig[16000];
 * // 1-second linear chirp from 200 Hz to 4000 Hz at 16 kHz sample rate
 * MD_chirp_linear(sig, 16000, 1.0, 200.0, 4000.0, 16000.0);
 * @endcode
 */
void MD_chirp_linear(double *output, unsigned N, double amplitude,
                     double f_start, double f_end, double sample_rate);

/**
 * Generate a logarithmic chirp (swept sine with exponentially increasing
 * frequency).
 *
 * The instantaneous frequency sweeps exponentially from @p f_start to
 * @p f_end over N samples:
 *
 *   f(t) = f_start * (f_end / f_start)^(t / T)
 *
 * where T = (N-1) / sample_rate is the sweep duration.  The output is:
 *
 *   output[i] = amplitude * sin(2*pi * f_start * T * (r^(t/T) - 1) / ln(r))
 *
 * where r = f_end / f_start.
 *
 * A logarithmic chirp spends equal time per octave, making it ideal for
 * measuring systems whose behaviour is best described on a log-frequency
 * axis (e.g. audio equaliser response).
 *
 * @param output      Output buffer of length N (caller-allocated).
 * @param N           Number of samples to generate.  Must be > 0.
 * @param amplitude   Peak amplitude of the chirp (e.g. 1.0).
 * @param f_start     Starting frequency in Hz.  Must be > 0.
 * @param f_end       Ending frequency in Hz.  Must be > 0 and != f_start.
 * @param sample_rate Sampling rate in Hz.  Must be > 0.
 *
 * @code
 * double sig[44100];
 * // 1-second log chirp from 20 Hz to 20 kHz at 44.1 kHz sample rate
 * MD_chirp_log(sig, 44100, 1.0, 20.0, 20000.0, 44100.0);
 * @endcode
 */
void MD_chirp_log(double *output, unsigned N, double amplitude,
                  double f_start, double f_end, double sample_rate);

/**
 * Generate a square wave.
 *
 * Fills the output buffer with a bipolar square wave that alternates
 * between +amplitude and −amplitude at the given frequency:
 *
 * \f[
 * x[n] = \begin{cases}
 *   +A  & 0 < \phi < \pi \\
 *   -A  & \pi < \phi < 2\pi \\
 *    0  & \phi = 0 \text{ or } \phi = \pi
 * \end{cases}
 * \f]
 *
 * where \f$\phi = 2\pi f n / f_s \pmod{2\pi}\f$.
 *
 * A square wave contains only odd harmonics (1f, 3f, 5f, …) whose
 * amplitudes decay as 1/k — a classic demonstration of the Fourier
 * series and the Gibbs phenomenon.
 *
 * @param output      Output buffer of length N (caller-allocated).
 * @param N           Number of samples to generate.  Must be > 0.
 * @param amplitude   Peak amplitude of the square wave (e.g. 1.0).
 * @param freq        Frequency in Hz.
 * @param sample_rate Sampling rate in Hz.  Must be > 0.
 *
 * @code
 * double sig[1024];
 * MD_square_wave(sig, 1024, 1.0, 440.0, 44100.0);
 * @endcode
 */
void MD_square_wave(double *output, unsigned N, double amplitude,
                    double freq, double sample_rate);

/**
 * Generate a sawtooth wave.
 *
 * Fills the output buffer with a sawtooth wave that ramps linearly
 * from −amplitude to +amplitude over each period:
 *
 * \f[
 * x[n] = A \left(\frac{\phi}{\pi} - 1\right)
 * \f]
 *
 * where \f$\phi = 2\pi f n / f_s \pmod{2\pi}\f$.
 *
 * A sawtooth wave contains all integer harmonics (1f, 2f, 3f, …)
 * whose amplitudes decay as 1/k — useful for demonstrating richer
 * harmonic content compared to the square wave's odd-only series.
 *
 * @param output      Output buffer of length N (caller-allocated).
 * @param N           Number of samples to generate.  Must be > 0.
 * @param amplitude   Peak amplitude of the sawtooth wave (e.g. 1.0).
 * @param freq        Frequency in Hz.
 * @param sample_rate Sampling rate in Hz.  Must be > 0.
 *
 * @code
 * double sig[1024];
 * MD_sawtooth_wave(sig, 1024, 1.0, 440.0, 44100.0);
 * @endcode
 */
void MD_sawtooth_wave(double *output, unsigned N, double amplitude,
                      double freq, double sample_rate);

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
 * @param ent         If non-null, receives the normalised entropy of the
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
