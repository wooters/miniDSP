/**
 * @file minidsp_generators.c
 * @brief Signal generators: sine, white noise, impulse, chirps, square,
 *        and sawtooth waves.
 * @author Chuck Wooters <wooters@hey.com>
 * @copyright 2013 International Computer Science Institute
 */

#include "minidsp.h"

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

void MD_white_noise(double *output, unsigned N, double amplitude,
                    unsigned seed)
{
    assert(output != nullptr);
    assert(N > 0);

    /* 64-bit LCG (Knuth MMIX constants) for uniform doubles in (0, 1).
     * Self-contained so we don't need POSIX drand48. */
    unsigned long long state = (unsigned long long)seed;

    /* Box-Muller transform: convert pairs of uniform random numbers
     * into pairs of Gaussian random numbers (mean 0, std dev 1),
     * then scale by the requested amplitude.
     *
     * For each pair (u1, u2) drawn from (0, 1]:
     *   z0 = sqrt(-2 ln u1) * cos(2 pi u2)
     *   z1 = sqrt(-2 ln u1) * sin(2 pi u2) */
    unsigned i = 0;
    while (i + 1 < N) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u1 = (double)(state >> 11) * 0x1.0p-53;
        if (u1 < DBL_MIN) u1 = DBL_MIN;
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u2 = (double)(state >> 11) * 0x1.0p-53;

        double r = sqrt(-2.0 * log(u1));
        double theta = 2.0 * M_PI * u2;
        output[i]     = amplitude * r * cos(theta);
        output[i + 1] = amplitude * r * sin(theta);
        i += 2;
    }
    /* If N is odd, generate one extra pair and use only the first value */
    if (i < N) {
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u1 = (double)(state >> 11) * 0x1.0p-53;
        if (u1 < DBL_MIN) u1 = DBL_MIN;
        state = state * 6364136223846793005ULL + 1442695040888963407ULL;
        double u2 = (double)(state >> 11) * 0x1.0p-53;
        output[i] = amplitude * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    }
}

void MD_impulse(double *output, unsigned N, double amplitude, unsigned position)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(position < N);
    memset(output, 0, N * sizeof(double));
    output[position] = amplitude;
}

void MD_chirp_linear(double *output, unsigned N, double amplitude,
                     double f_start, double f_end, double sample_rate)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);

    if (N == 1) {
        output[0] = 0.0;
        return;
    }

    double T = (double)(N - 1) / sample_rate;
    double chirp_rate = (f_end - f_start) / T;

    for (unsigned i = 0; i < N; i++) {
        double t = (double)i / sample_rate;
        double phase = 2.0 * M_PI * (f_start * t + 0.5 * chirp_rate * t * t);
        output[i] = amplitude * sin(phase);
    }
}

void MD_chirp_log(double *output, unsigned N, double amplitude,
                  double f_start, double f_end, double sample_rate)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);
    assert(f_start > 0.0);
    assert(f_end > 0.0);
    assert(f_start != f_end);

    if (N == 1) {
        output[0] = 0.0;
        return;
    }

    double T = (double)(N - 1) / sample_rate;
    double ratio = f_end / f_start;
    double log_ratio = log(ratio);

    for (unsigned i = 0; i < N; i++) {
        double t = (double)i / sample_rate;
        double phase = 2.0 * M_PI * f_start * T
                     * (pow(ratio, t / T) - 1.0) / log_ratio;
        output[i] = amplitude * sin(phase);
    }
}

void MD_square_wave(double *output, unsigned N, double amplitude,
                    double freq, double sample_rate)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);
    double phase_step = 2.0 * M_PI * freq / sample_rate;
    for (unsigned i = 0; i < N; i++) {
        double phase = fmod(phase_step * (double)i, 2.0 * M_PI);
        if (phase < 0.0) phase += 2.0 * M_PI;
        if (phase == 0.0 || phase == M_PI)
            output[i] = 0.0;
        else if (phase < M_PI)
            output[i] = amplitude;
        else
            output[i] = -amplitude;
    }
}

void MD_sawtooth_wave(double *output, unsigned N, double amplitude,
                      double freq, double sample_rate)
{
    assert(output != nullptr);
    assert(N > 0);
    assert(sample_rate > 0.0);
    double phase_step = 2.0 * M_PI * freq / sample_rate;
    for (unsigned i = 0; i < N; i++) {
        double phase = fmod(phase_step * (double)i, 2.0 * M_PI);
        if (phase < 0.0) phase += 2.0 * M_PI;
        output[i] = amplitude * (phase / M_PI - 1.0);
    }
}
