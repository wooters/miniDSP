/**
 * @file minidsp_resample.c
 * @brief Polyphase sinc resampler for offline sample rate conversion.
 */

#include "minidsp.h"
#include "minidsp_internal.h"

/** Number of sub-phases in the polyphase filter table. */
#define RESAMPLE_NUM_PHASES 512

unsigned MD_resample_output_len(unsigned input_len,
                                double in_rate, double out_rate)
{
    MD_CHECK(input_len > 0, MD_ERR_INVALID_SIZE, "input_len must be > 0", 0);
    MD_CHECK(in_rate > 0.0, MD_ERR_INVALID_RANGE, "in_rate must be > 0", 0);
    MD_CHECK(out_rate > 0.0, MD_ERR_INVALID_RANGE, "out_rate must be > 0", 0);
    return (unsigned)ceil((double)input_len * out_rate / in_rate);
}

unsigned MD_resample(const double *input, unsigned input_len,
                     double *output, unsigned max_output_len,
                     double in_rate, double out_rate,
                     unsigned num_zero_crossings, double kaiser_beta)
{
    MD_CHECK(input != NULL, MD_ERR_NULL_POINTER, "input must not be NULL", 0);
    MD_CHECK(output != NULL, MD_ERR_NULL_POINTER, "output must not be NULL", 0);
    MD_CHECK(input_len > 0, MD_ERR_INVALID_SIZE, "input_len must be > 0", 0);
    MD_CHECK(in_rate > 0.0, MD_ERR_INVALID_RANGE, "in_rate must be > 0", 0);
    MD_CHECK(out_rate > 0.0, MD_ERR_INVALID_RANGE, "out_rate must be > 0", 0);
    MD_CHECK(num_zero_crossings > 0, MD_ERR_INVALID_SIZE, "num_zero_crossings must be > 0", 0);

    unsigned expected_len = MD_resample_output_len(input_len, in_rate, out_rate);
    MD_CHECK(max_output_len >= expected_len, MD_ERR_INVALID_SIZE,
             "max_output_len too small for expected output", 0);

    unsigned num_phases = RESAMPLE_NUM_PHASES;
    unsigned taps_per_phase = 2 * num_zero_crossings;
    unsigned table_size = (num_phases + 1) * taps_per_phase;

    /* Build polyphase filter table.
     * table[phase][tap] contains the windowed-sinc value for that
     * sub-phase offset and tap index. We build num_phases+1 phases
     * so linear interpolation between phase[p] and phase[p+1] works
     * at the boundary. */
    double *table = malloc(table_size * sizeof(double));
    MD_CHECK(table != NULL, MD_ERR_ALLOC_FAILED, "malloc failed", 0);

    /* Anti-aliasing: scale cutoff for downsampling */
    double ratio = out_rate / in_rate;
    double cutoff_ratio = (ratio < 1.0) ? ratio : 1.0;

    for (unsigned p = 0; p <= num_phases; p++) {
        double frac = (double)p / (double)num_phases;
        double *phase_row = table + p * taps_per_phase;

        for (unsigned t = 0; t < taps_per_phase; t++) {
            int tap_offset = (int)t - (int)num_zero_crossings;
            double x = (double)tap_offset + frac;

            /* Windowed sinc: sinc(cutoff * x) * kaiser(x / half_width) */
            double sinc_val = MD_sinc(cutoff_ratio * x);

            /* Kaiser window over the full filter span */
            double half_width = (double)num_zero_crossings;
            double w;
            if (fabs(x) >= half_width) {
                w = 0.0;
            } else {
                double r = x / half_width;
                double arg = 1.0 - r * r;
                if (arg < 0.0) arg = 0.0;
                w = MD_bessel_i0(kaiser_beta * sqrt(arg))
                  / MD_bessel_i0(kaiser_beta);
            }

            phase_row[t] = cutoff_ratio * sinc_val * w;
        }
    }

    /* Resample: for each output sample, compute fractional input position,
     * select sub-phase with linear interpolation, dot with input. */
    unsigned out_len = expected_len;

    for (unsigned n = 0; n < out_len; n++) {
        double in_pos = (double)n / ratio;
        int in_idx = (int)floor(in_pos);
        double frac = in_pos - (double)in_idx;

        /* Map fractional position to phase index */
        double phase_exact = frac * (double)num_phases;
        unsigned phase_lo = (unsigned)phase_exact;
        double alpha = phase_exact - (double)phase_lo;

        if (phase_lo >= num_phases) {
            phase_lo = num_phases - 1;
            alpha = 1.0;
        }

        const double *row_lo = table + phase_lo * taps_per_phase;
        const double *row_hi = table + (phase_lo + 1) * taps_per_phase;

        double acc = 0.0;
        for (unsigned t = 0; t < taps_per_phase; t++) {
            int si = in_idx + (int)t - (int)num_zero_crossings;

            /* Clamp to input boundaries */
            double sample;
            if (si < 0 || (unsigned)si >= input_len) {
                sample = 0.0;
            } else {
                sample = input[si];
            }

            double coeff = row_lo[t] + alpha * (row_hi[t] - row_lo[t]);
            acc += sample * coeff;
        }

        output[n] = acc;
    }

    free(table);
    return out_len;
}
