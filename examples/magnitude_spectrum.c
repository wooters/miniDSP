/**
 * @file magnitude_spectrum.c
 * @brief Example: compute and export the magnitude spectrum of a test signal.
 *
 * This program demonstrates MD_magnitude_spectrum() by:
 *   1. Generating a signal with three sinusoidal components (440 Hz,
 *      1000 Hz, and 2500 Hz) plus a DC offset.
 *   2. Applying a Hanning window to reduce spectral leakage.
 *   3. Computing the magnitude spectrum via MD_magnitude_spectrum().
 *   4. Writing the results to a CSV file for plotting.
 *
 * Build (from the repository root):
 *   gcc-14 -std=c23 -Wall -Wextra -pedantic -O2 -I.. \
 *       examples/magnitude_spectrum.c -L. -lminidsp -lfftw3 -lm \
 *       -o examples/magnitude_spectrum
 *
 * Run:
 *   ./examples/magnitude_spectrum
 *   python3 examples/plot_spectrum.py   # generates magnitude_spectrum.png
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"

int main(void)
{
    /* ------------------------------------------------------------------
     * Signal parameters
     * ----------------------------------------------------------------*/
    const unsigned   N           = 4096;        /* FFT size (samples)       */
    const double     sample_rate = 16000.0;     /* Hz                       */
    const unsigned   num_bins    = N / 2 + 1;   /* unique frequency bins    */

    /* Three tones at different amplitudes, plus a small DC offset */
    const double     freq1 =  440.0;   /* A4 note               */
    const double     amp1  =    1.0;
    const double     freq2 = 1000.0;   /* 1 kHz reference tone  */
    const double     amp2  =    0.6;
    const double     freq3 = 2500.0;   /* higher partial        */
    const double     amp3  =    0.3;
    const double     dc    =    0.1;   /* small DC offset       */

    /* ------------------------------------------------------------------
     * Generate the test signal
     * ----------------------------------------------------------------*/
    double *signal  = malloc(N * sizeof(double));
    double *windowed = malloc(N * sizeof(double));
    double *window  = malloc(N * sizeof(double));
    double *mag     = malloc(num_bins * sizeof(double));

    if (!signal || !windowed || !window || !mag) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    for (unsigned i = 0; i < N; i++) {
        double t = (double)i / sample_rate;
        signal[i] = dc
                  + amp1 * sin(2.0 * M_PI * freq1 * t)
                  + amp2 * sin(2.0 * M_PI * freq2 * t)
                  + amp3 * sin(2.0 * M_PI * freq3 * t);
    }

    /* ------------------------------------------------------------------
     * Apply a Hanning window to reduce spectral leakage
     * ----------------------------------------------------------------*/
    MD_Gen_Hann_Win(window, N);
    for (unsigned i = 0; i < N; i++) {
        windowed[i] = signal[i] * window[i];
    }

    /* ------------------------------------------------------------------
     * Compute the magnitude spectrum
     * ----------------------------------------------------------------*/
    MD_magnitude_spectrum(windowed, N, mag);

    /* Convert to single-sided amplitude spectrum:
     *   - Divide all bins by N
     *   - Multiply interior bins by 2 (because negative frequencies
     *     are folded into the positive side)
     *   - DC and Nyquist bins appear only once, so no doubling */
    for (unsigned k = 0; k < num_bins; k++) {
        mag[k] /= (double)N;
        if (k > 0 && k < N / 2) {
            mag[k] *= 2.0;
        }
    }

    /* ------------------------------------------------------------------
     * Write results to CSV
     * ----------------------------------------------------------------*/
    const char *csv_file = "examples/magnitude_spectrum.csv";
    FILE *fp = fopen(csv_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", csv_file);
        return 1;
    }

    fprintf(fp, "bin,frequency_hz,magnitude\n");
    for (unsigned k = 0; k < num_bins; k++) {
        double freq_hz = (double)k * sample_rate / (double)N;
        fprintf(fp, "%u,%.4f,%.8f\n", k, freq_hz, mag[k]);
    }
    fclose(fp);

    printf("Wrote %u frequency bins to %s\n", num_bins, csv_file);
    printf("Signal: %.0f Hz + %.0f Hz + %.0f Hz + DC offset\n",
           freq1, freq2, freq3);
    printf("Sample rate: %.0f Hz, FFT size: %u\n", sample_rate, N);

    /* ------------------------------------------------------------------
     * Cleanup
     * ----------------------------------------------------------------*/
    free(mag);
    free(window);
    free(windowed);
    free(signal);
    MD_shutdown();

    return 0;
}
