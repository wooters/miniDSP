/**
 * @file minidsp_steg.c
 * @brief Audio steganography: hide secret messages inside audio signals.
 *
 * Two methods are provided:
 *
 *   - **LSB (Least Significant Bit)** — encodes message bits in the lowest
 *     bit of a 16-bit PCM representation of each sample.  High capacity,
 *     negligible audible distortion, but fragile to any signal processing.
 *
 *   - **Frequency-band modulation** — encodes bits as the presence or
 *     absence of a near-ultrasonic tone (default 18.5 kHz / 19.5 kHz BFSK)
 *     in short time chips.  Lower capacity but survives mild resampling
 *     and compression.  Requires sample_rate >= 44100 Hz so the carrier
 *     frequencies remain below Nyquist.
 *
 * Both methods prepend a 32-bit little-endian length header (message byte
 * count) before the payload, enabling blind decoding without knowing the
 * original message length.
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Shared constants
 * ----------------------------------------------------------------------- */

/** Bits needed for the message-length header (32-bit unsigned). */
#define HEADER_BITS 32

/* -----------------------------------------------------------------------
 * LSB steganography
 *
 * The host signal is in [-1, 1] double.  We map it to a 16-bit PCM
 * integer, flip the LSB to carry one message bit per sample, then map
 * back.  The distortion introduced is at most ±1/32768 ≈ −90 dB.
 *
 * Bit layout:  [32 bits: length (LE)] [8*length bits: message bytes]
 * Each bit uses one sample.
 * ----------------------------------------------------------------------- */

/** Clamp a double to [-1, 1]. */
static double clamp(double x)
{
    if (x >  1.0) return  1.0;
    if (x < -1.0) return -1.0;
    return x;
}

/** Convert a double sample in [-1, 1] to a signed 16-bit PCM value.
 *  Uses rounding (not truncation) so the conversion survives a
 *  double → float → double WAV roundtrip. */
static int double_to_pcm16(double sample)
{
    double s = clamp(sample);
    int pcm = (s >= 0.0) ? (int)(s * 32767.0 + 0.5)
                          : (int)(s * 32767.0 - 0.5);
    if (pcm >  32767) pcm =  32767;
    if (pcm < -32767) pcm = -32767;
    return pcm;
}

/** Convert a signed 16-bit PCM value back to double in [-1, 1].
 *  The result is snapped to float precision so that a WAV write/read
 *  cycle (which stores IEEE float) preserves the exact value. */
static double pcm16_to_double(int pcm)
{
    return (double)(float)((double)pcm / 32767.0);
}

static unsigned lsb_capacity(unsigned signal_len)
{
    if (signal_len <= HEADER_BITS)
        return 0;
    return (signal_len - HEADER_BITS) / 8;
}

static unsigned lsb_encode(const double *host, double *output,
                           unsigned signal_len, const char *message)
{
    unsigned msg_len = (unsigned)strlen(message);
    unsigned capacity = lsb_capacity(signal_len);
    if (msg_len == 0 || capacity == 0)
        return 0;
    if (msg_len > capacity)
        msg_len = capacity;

    unsigned total_bits = HEADER_BITS + msg_len * 8;

    /* Copy host to output. */
    memcpy(output, host, signal_len * sizeof(double));

    /* Encode the 32-bit length header (little-endian). */
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = (msg_len >> i) & 1;
        int pcm = double_to_pcm16(output[i]);
        /* Clear the LSB and set it to our message bit. */
        pcm = (pcm & ~1) | (int)bit;
        output[i] = pcm16_to_double(pcm);
    }

    /* Encode the message bytes. */
    for (unsigned b = 0; b < msg_len; b++) {
        unsigned char ch = (unsigned char)message[b];
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned sample_idx = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = (ch >> bit_idx) & 1;
            int pcm = double_to_pcm16(output[sample_idx]);
            pcm = (pcm & ~1) | (int)bit;
            output[sample_idx] = pcm16_to_double(pcm);
        }
    }

    (void)total_bits;
    return msg_len;
}

static unsigned lsb_decode(const double *stego, unsigned signal_len,
                           char *message_out, unsigned max_msg_len)
{
    if (signal_len <= HEADER_BITS || max_msg_len == 0)
        return 0;

    /* Read the 32-bit length header. */
    unsigned msg_len = 0;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        int pcm = double_to_pcm16(stego[i]);
        unsigned bit = (unsigned)(pcm & 1);
        msg_len |= (bit << i);
    }

    /* Sanity check. */
    unsigned capacity = lsb_capacity(signal_len);
    if (msg_len == 0 || msg_len > capacity)
        return 0;

    unsigned decode_len = msg_len;
    if (decode_len >= max_msg_len)
        decode_len = max_msg_len - 1;

    for (unsigned b = 0; b < decode_len; b++) {
        unsigned char ch = 0;
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned sample_idx = HEADER_BITS + b * 8 + bit_idx;
            int pcm = double_to_pcm16(stego[sample_idx]);
            unsigned bit = (unsigned)(pcm & 1);
            ch |= (unsigned char)(bit << bit_idx);
        }
        message_out[b] = (char)ch;
    }
    message_out[decode_len] = '\0';
    return decode_len;
}

/* -----------------------------------------------------------------------
 * Frequency-band steganography (BFSK)
 *
 * Uses Binary Frequency-Shift Keying in the near-ultrasonic range:
 *   - bit 0 → tone at FREQ_LO (18.5 kHz)
 *   - bit 1 → tone at FREQ_HI (19.5 kHz)
 *
 * Each bit occupies CHIP_MS milliseconds of audio.  The tone is mixed
 * additively at a low amplitude so it's inaudible to most listeners.
 *
 * Decoding correlates each chip against both carriers and picks the
 * stronger one.
 *
 * Bit layout: [32 bits: length (LE)] [8*length bits: message bytes]
 * ----------------------------------------------------------------------- */

#define FREQ_LO    18500.0   /**< Carrier for bit 0 (Hz). */
#define FREQ_HI    19500.0   /**< Carrier for bit 1 (Hz). */
#define CHIP_MS    3.0       /**< Duration of one bit chip (ms). */
#define TONE_AMP   0.02      /**< Additive tone amplitude. */

static unsigned chip_samples(double sample_rate)
{
    return (unsigned)(CHIP_MS * sample_rate / 1000.0);
}

static unsigned freq_capacity(unsigned signal_len, double sample_rate)
{
    unsigned cs = chip_samples(sample_rate);
    if (cs == 0) return 0;
    unsigned total_chips = signal_len / cs;
    if (total_chips <= HEADER_BITS)
        return 0;
    return (total_chips - HEADER_BITS) / 8;
}

static unsigned freq_encode(const double *host, double *output,
                            unsigned signal_len, double sample_rate,
                            const char *message)
{
    unsigned msg_len = (unsigned)strlen(message);
    unsigned capacity = freq_capacity(signal_len, sample_rate);
    if (msg_len == 0 || capacity == 0)
        return 0;
    if (msg_len > capacity)
        msg_len = capacity;

    unsigned cs = chip_samples(sample_rate);

    /* Copy host to output. */
    memcpy(output, host, signal_len * sizeof(double));

    /* Helper: add a tone burst for one bit at the given chip index. */
    #define ENCODE_BIT(chip_idx, bit_val) do {                      \
        double freq = (bit_val) ? FREQ_HI : FREQ_LO;               \
        unsigned start = (chip_idx) * cs;                           \
        for (unsigned s = 0; s < cs && (start + s) < signal_len; s++) { \
            double t = (double)s / sample_rate;                     \
            output[start + s] += TONE_AMP * sin(2.0 * M_PI * freq * t); \
        }                                                           \
    } while (0)

    /* Encode 32-bit length header. */
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = (msg_len >> i) & 1;
        ENCODE_BIT(i, bit);
    }

    /* Encode message bytes. */
    for (unsigned b = 0; b < msg_len; b++) {
        unsigned char ch = (unsigned char)message[b];
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned chip = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = (ch >> bit_idx) & 1;
            ENCODE_BIT(chip, bit);
        }
    }

    #undef ENCODE_BIT
    return msg_len;
}

/** Decode one bit by correlating a chip against both BFSK carriers. */
static unsigned decode_one_bit(const double *stego, unsigned signal_len,
                               double sample_rate, unsigned chip_idx,
                               unsigned cs)
{
    unsigned start = chip_idx * cs;
    double corr_lo = 0.0, corr_hi = 0.0;
    for (unsigned s = 0; s < cs && (start + s) < signal_len; s++) {
        double t = (double)s / sample_rate;
        corr_lo += stego[start + s] * sin(2.0 * M_PI * FREQ_LO * t);
        corr_hi += stego[start + s] * sin(2.0 * M_PI * FREQ_HI * t);
    }
    return (fabs(corr_hi) > fabs(corr_lo)) ? 1u : 0u;
}

static unsigned freq_decode(const double *stego, unsigned signal_len,
                            double sample_rate,
                            char *message_out, unsigned max_msg_len)
{
    unsigned cs = chip_samples(sample_rate);
    if (cs == 0 || max_msg_len == 0)
        return 0;

    unsigned total_chips = signal_len / cs;
    if (total_chips <= HEADER_BITS)
        return 0;

    /* Read the 32-bit length header. */
    unsigned msg_len = 0;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = decode_one_bit(stego, signal_len, sample_rate, i, cs);
        msg_len |= (bit << i);
    }

    /* Sanity check. */
    unsigned capacity = freq_capacity(signal_len, sample_rate);
    if (msg_len == 0 || msg_len > capacity)
        return 0;

    unsigned decode_len = msg_len;
    if (decode_len >= max_msg_len)
        decode_len = max_msg_len - 1;

    for (unsigned b = 0; b < decode_len; b++) {
        unsigned char ch = 0;
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned chip = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = decode_one_bit(stego, signal_len, sample_rate,
                                          chip, cs);
            ch |= (unsigned char)(bit << bit_idx);
        }
        message_out[b] = (char)ch;
    }
    message_out[decode_len] = '\0';
    return decode_len;
}

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

unsigned MD_steg_capacity(unsigned signal_len, double sample_rate, int method)
{
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND);

    if (method == MD_STEG_LSB)
        return lsb_capacity(signal_len);
    else
        return freq_capacity(signal_len, sample_rate);
}

unsigned MD_steg_encode(const double *host, double *output,
                        unsigned signal_len, double sample_rate,
                        const char *message, int method)
{
    assert(host != nullptr);
    assert(output != nullptr);
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(message != nullptr);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND);

    if (method == MD_STEG_FREQ_BAND)
        assert(sample_rate >= 40000.0 &&
               "frequency-band steganography requires sample_rate >= 40 kHz");

    if (method == MD_STEG_LSB)
        return lsb_encode(host, output, signal_len, message);
    else
        return freq_encode(host, output, signal_len, sample_rate, message);
}

unsigned MD_steg_decode(const double *stego, unsigned signal_len,
                        double sample_rate,
                        char *message_out, unsigned max_msg_len,
                        int method)
{
    assert(stego != nullptr);
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(message_out != nullptr);
    assert(max_msg_len > 0);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND);

    if (method == MD_STEG_LSB)
        return lsb_decode(stego, signal_len, message_out, max_msg_len);
    else
        return freq_decode(stego, signal_len, sample_rate,
                           message_out, max_msg_len);
}
