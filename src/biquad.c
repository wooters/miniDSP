/**
 * @file biquad.c
 * @brief Biquad (second-order IIR) filter implementation.
 *
 * Original implementation by Tom St Denis, based on:
 *   "Cookbook formulae for audio EQ biquad filter coefficients"
 *   by Robert Bristow-Johnson  (pbjrbj@viconet.com)
 *
 * Reference:
 *   https://www.w3.org/2011/audio/audio-eq-cookbook.html
 *
 * This code is in the public domain.  See the original header:
 *   http://www.musicdsp.org/files/biquad.c
 *
 * -----------------------------------------------------------------------
 *
 * How a biquad filter works (for students):
 *
 * A biquad filter is defined by the "difference equation":
 *
 *   y[n] = (b0/a0)*x[n] + (b1/a0)*x[n-1] + (b2/a0)*x[n-2]
 *                        - (a1/a0)*y[n-1] - (a2/a0)*y[n-2]
 *
 * where x[n] is the current input sample, y[n] is the current output,
 * and the subscripted values are previous samples.  The six coefficients
 * (a0, a1, a2, b0, b1, b2) completely determine the filter's behaviour.
 *
 * Different choices of coefficients produce different filter types
 * (low-pass, high-pass, etc.).  The BiQuad_new() function computes
 * the right coefficients from user-friendly parameters like frequency,
 * bandwidth, and gain.
 */
#include "biquad.h"

/**
 * Process one input sample through the biquad filter.
 *
 * This is the "hot loop" function -- it gets called once for every
 * single audio sample (e.g. 44,100 times per second at CD quality).
 * That is why it must be very simple and fast.
 *
 * The five multiply-accumulate operations implement:
 *   result = a0*x[n] + a1*x[n-1] + a2*x[n-2] - a3*y[n-1] - a4*y[n-2]
 *
 * After computing the result, we shift the delay line: the current
 * sample becomes x[n-1] for the next call, and x[n-1] becomes x[n-2].
 * Same for the output side.
 */
smp_type BiQuad(smp_type sample, biquad *b)
{
    smp_type result;

    /* Apply the difference equation */
    result = b->a0 * sample
           + b->a1 * b->x1
           + b->a2 * b->x2
           - b->a3 * b->y1
           - b->a4 * b->y2;

    /* Shift the input delay line */
    b->x2 = b->x1;
    b->x1 = sample;

    /* Shift the output delay line */
    b->y2 = b->y1;
    b->y1 = result;

    return result;
}

/**
 * Create a new biquad filter with the specified characteristics.
 *
 * The function computes the six transfer-function coefficients
 * (b0, b1, b2, a0, a1, a2) from the desired filter type and
 * parameters, then normalises them by dividing through by a0.
 *
 * @param type       Filter type (LPF, HPF, BPF, NOTCH, PEQ, LSH, HSH).
 * @param dbGain     Gain in dB.  Only matters for PEQ, LSH, and HSH.
 * @param freq       Centre/corner frequency in Hz.
 * @param srate      Sampling rate in Hz.
 * @param bandwidth  Bandwidth in octaves.
 * @return           Pointer to a new biquad struct, or NULL on error.
 *
 * Note: the caller is responsible for calling free() on the returned
 * pointer when the filter is no longer needed.
 */
biquad *BiQuad_new(int type, smp_type dbGain,
                   smp_type freq, smp_type srate,
                   smp_type bandwidth)
{
    biquad *b;
    smp_type A, omega, sn, cs, alpha, beta;
    smp_type a0, a1, a2, b0, b1, b2;

    b = malloc(sizeof(biquad));
    if (b == NULL)
        return NULL;

    /* --- Compute intermediate variables ---
     *
     * A     = linear amplitude from dB gain (used by PEQ and shelf filters)
     * omega = angular frequency of the centre/corner in radians per sample
     * sn    = sin(omega)
     * cs    = cos(omega)
     * alpha = controls the bandwidth (Q factor)
     * beta  = helper for shelf filter coefficients
     */
    A     = pow(10.0, dbGain / 40.0);
    omega = 2.0 * M_PI * freq / srate;
    sn    = sin(omega);
    cs    = cos(omega);
    alpha = sn * sinh(M_LN2 / 2.0 * bandwidth * omega / sn);
    /* beta * sn substitutes for the cookbook's 2*sqrt(A)*alpha in shelf
     * filters (LSH/HSH).  This is mathematically equivalent to fixing the
     * shelf slope parameter S = 1 (steepest slope without resonance).
     * Because of this, the 'bandwidth' parameter has no effect on LSH/HSH. */
    beta  = sqrt(A + A);  /* = sqrt(2*A) */

    switch (type) {

    /* --- Low-Pass Filter ---
     * Passes frequencies below 'freq', attenuates higher ones.
     * Think of it as removing treble from an audio signal. */
    case LPF:
        b0 = (1.0 - cs) / 2.0;
        b1 =  1.0 - cs;
        b2 = (1.0 - cs) / 2.0;
        a0 =  1.0 + alpha;
        a1 = -2.0 * cs;
        a2 =  1.0 - alpha;
        break;

    /* --- High-Pass Filter ---
     * Passes frequencies above 'freq', attenuates lower ones.
     * Useful for removing low-frequency rumble (wind noise, etc.). */
    case HPF:
        b0 = (1.0 + cs) / 2.0;
        b1 = -(1.0 + cs);
        b2 = (1.0 + cs) / 2.0;
        a0 =  1.0 + alpha;
        a1 = -2.0 * cs;
        a2 =  1.0 - alpha;
        break;

    /* --- Band-Pass Filter ---
     * Passes a range of frequencies centred on 'freq'.
     * Width is controlled by 'bandwidth' (in octaves). */
    case BPF:
        b0 =  alpha;
        b1 =  0.0;
        b2 = -alpha;
        a0 =  1.0 + alpha;
        a1 = -2.0 * cs;
        a2 =  1.0 - alpha;
        break;

    /* --- Notch (Band-Reject) Filter ---
     * Removes a narrow band of frequencies around 'freq'.
     * Commonly used to eliminate 50/60 Hz mains hum. */
    case NOTCH:
        b0 =  1.0;
        b1 = -2.0 * cs;
        b2 =  1.0;
        a0 =  1.0 + alpha;
        a1 = -2.0 * cs;
        a2 =  1.0 - alpha;
        break;

    /* --- Peaking EQ Filter ---
     * Boosts or cuts a band of frequencies by 'dbGain' dB.
     * This is the "parametric EQ" knob on a mixing console. */
    case PEQ:
        b0 =  1.0 + (alpha * A);
        b1 = -2.0 * cs;
        b2 =  1.0 - (alpha * A);
        a0 =  1.0 + (alpha / A);
        a1 = -2.0 * cs;
        a2 =  1.0 - (alpha / A);
        break;

    /* --- Low Shelf Filter ---
     * Boosts or cuts all frequencies below 'freq' by 'dbGain' dB.
     * Like a bass knob on a stereo. */
    case LSH:
        b0 = A * ((A + 1.0) - (A - 1.0) * cs + beta * sn);
        b1 = 2.0 * A * ((A - 1.0) - (A + 1.0) * cs);
        b2 = A * ((A + 1.0) - (A - 1.0) * cs - beta * sn);
        a0 = (A + 1.0) + (A - 1.0) * cs + beta * sn;
        a1 = -2.0 * ((A - 1.0) + (A + 1.0) * cs);
        a2 = (A + 1.0) + (A - 1.0) * cs - beta * sn;
        break;

    /* --- High Shelf Filter ---
     * Boosts or cuts all frequencies above 'freq' by 'dbGain' dB.
     * Like a treble knob on a stereo. */
    case HSH:
        b0 = A * ((A + 1.0) + (A - 1.0) * cs + beta * sn);
        b1 = -2.0 * A * ((A - 1.0) + (A + 1.0) * cs);
        b2 = A * ((A + 1.0) + (A - 1.0) * cs - beta * sn);
        a0 = (A + 1.0) - (A - 1.0) * cs + beta * sn;
        a1 = 2.0 * ((A - 1.0) - (A + 1.0) * cs);
        a2 = (A + 1.0) - (A - 1.0) * cs - beta * sn;
        break;

    default:
        free(b);
        return NULL;
    }

    /* Normalise all coefficients by dividing through by a0.
     * This means the filter equation becomes:
     *   y[n] = a0*x[n] + a1*x[n-1] + a2*x[n-2] - a3*y[n-1] - a4*y[n-2]
     * where a0..a4 are the normalised coefficients stored in the struct. */
    b->a0 = b0 / a0;
    b->a1 = b1 / a0;
    b->a2 = b2 / a0;
    b->a3 = a1 / a0;
    b->a4 = a2 / a0;

    /* Clear the delay lines (no previous samples yet) */
    b->x1 = b->x2 = 0.0;
    b->y1 = b->y2 = 0.0;

    return b;
}
