/**
 * @file resampler.c
 * @brief Example: resample a signal between common audio sample rates.
 *
 * Demonstrates:
 *   1. MD_resample_output_len() -- compute output buffer size
 *   2. MD_resample() -- polyphase sinc resampling
 *
 * Build and run (from repo root):
 *   make -C examples resampler
 *   ./examples/resampler
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"

int main(void)
{
    /* Generate a 440 Hz sine wave at 44100 Hz (1 second). */
    const unsigned N_in = 44100;
    const double in_rate = 44100.0;
    double *input = malloc(N_in * sizeof(double));
    MD_sine_wave(input, N_in, 1.0, 440.0, in_rate);

    /* Resample to 48000 Hz. */
    //! [resample-basic]
    const double out_rate = 48000.0;
    unsigned N_out = MD_resample_output_len(N_in, in_rate, out_rate);
    double *output = malloc(N_out * sizeof(double));

    unsigned n_written = MD_resample(input, N_in, output, N_out,
                                     in_rate, out_rate, 32, 10.0);
    //! [resample-basic]

    printf("Resampled %u samples at %.0f Hz -> %u samples at %.0f Hz\n",
           N_in, in_rate, n_written, out_rate);

    /* Verify: measure F0 of the output. */
    double f0 = MD_f0_autocorrelation(output, n_written, out_rate,
                                       100.0, 2000.0);
    printf("Input frequency: 440 Hz\n");
    printf("Output F0:       %.1f Hz\n", f0);

    free(output);
    free(input);
    return 0;
}
