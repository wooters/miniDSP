/**
 * @file minidsp_steg.c
 * @brief Audio steganography: hide secret messages inside audio signals.
 *
 * Three methods are provided:
 *
 *   - **LSB (Least Significant Bit)** — encodes message bits in the lowest
 *     bit of a 16-bit PCM representation of each sample.  High capacity,
 *     negligible audible distortion, but fragile to any signal processing.
 *
 *   - **Frequency-band modulation** — encodes bits as the presence or
 *     absence of a near-ultrasonic tone (default 18.5 kHz / 19.5 kHz BFSK)
 *     in short time chips.  Lower capacity but survives mild resampling
 *     and compression.  Requires sample_rate >= 40000 Hz so the carrier
 *     frequencies remain below Nyquist.
 *
 *   - **Spectrogram text (spectext)** — hybrid: LSB for machine-readable
 *     data + spectrogram text art in the 18-24 kHz ultrasonic band for
 *     visual verification.  Auto-upsamples to 48 kHz.  The visual channel
 *     uses a fixed 30 ms column width (240 ms per character).
 *
 * Both methods prepend a 32-bit little-endian header before the payload.
 * Bits 0–30 hold the message byte count; bit 31 is a payload type flag
 * (0 = text, 1 = binary).  This enables blind decoding without knowing the
 * original message length, and allows MD_steg_detect() to identify the
 * payload type.
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Shared constants
 * ----------------------------------------------------------------------- */

/** Bits needed for the message-length header (32-bit unsigned). */
#define HEADER_BITS 32

/** Bit 31 of the header: payload type flag (0 = text, 1 = binary). */
#define PAYLOAD_TYPE_BIT (1u << 31)

/** Mask to extract the message length from the raw header (bits 0–30). */
#define LENGTH_MASK 0x7FFFFFFFu

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
                           unsigned signal_len,
                           const unsigned char *data, unsigned data_len,
                           unsigned flags)
{
    unsigned capacity = lsb_capacity(signal_len);
    if (data_len == 0 || capacity == 0)
        return 0;
    if (data_len > capacity)
        data_len = capacity;

    /* Copy host to output. */
    memcpy(output, host, signal_len * sizeof(double));

    /* Encode the 32-bit length header (little-endian).
     * Bit 31 carries the payload type flag. */
    unsigned header = data_len | flags;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = (header >> i) & 1;
        int pcm = double_to_pcm16(output[i]);
        /* Clear the LSB and set it to our message bit. */
        pcm = (pcm & ~1) | (int)bit;
        output[i] = pcm16_to_double(pcm);
    }

    /* Encode the data bytes. */
    for (unsigned b = 0; b < data_len; b++) {
        unsigned char ch = data[b];
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned sample_idx = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = (ch >> bit_idx) & 1;
            int pcm = double_to_pcm16(output[sample_idx]);
            pcm = (pcm & ~1) | (int)bit;
            output[sample_idx] = pcm16_to_double(pcm);
        }
    }

    return data_len;
}

/** Read the raw 32-bit LSB header from the first HEADER_BITS samples. */
static unsigned lsb_read_header(const double *signal)
{
    unsigned raw = 0;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        int pcm = double_to_pcm16(signal[i]);
        unsigned bit = (unsigned)(pcm & 1);
        raw |= (bit << i);
    }
    return raw;
}

static unsigned lsb_decode(const double *stego, unsigned signal_len,
                           unsigned char *data_out, unsigned max_len)
{
    if (signal_len <= HEADER_BITS || max_len == 0)
        return 0;

    unsigned msg_len = lsb_read_header(stego) & LENGTH_MASK;

    /* Sanity check. */
    unsigned capacity = lsb_capacity(signal_len);
    if (msg_len == 0 || msg_len > capacity)
        return 0;

    unsigned decode_len = msg_len;
    if (decode_len > max_len)
        decode_len = max_len;

    for (unsigned b = 0; b < decode_len; b++) {
        unsigned char ch = 0;
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned sample_idx = HEADER_BITS + b * 8 + bit_idx;
            int pcm = double_to_pcm16(stego[sample_idx]);
            unsigned bit = (unsigned)(pcm & 1);
            ch |= (unsigned char)(bit << bit_idx);
        }
        data_out[b] = ch;
    }
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

static unsigned freq_capacity_cs(unsigned signal_len, unsigned cs)
{
    if (cs == 0) return 0;
    unsigned total_chips = signal_len / cs;
    if (total_chips <= HEADER_BITS)
        return 0;
    return (total_chips - HEADER_BITS) / 8;
}

static unsigned freq_capacity(unsigned signal_len, double sample_rate)
{
    return freq_capacity_cs(signal_len, chip_samples(sample_rate));
}

/* -----------------------------------------------------------------------
 * Cached BFSK carrier sine tables
 *
 * The sine lookup tables depend only on sample_rate (which determines
 * chip_samples).  We cache them across calls and only recompute when
 * the sample rate changes.  MD_shutdown() frees them via
 * md_steg_teardown() (declared in minidsp_internal.h).
 * ----------------------------------------------------------------------- */

static double  *_bfsk_sin_lo  = NULL;
static double  *_bfsk_sin_hi  = NULL;
static unsigned  _bfsk_cs     = 0;

static void _bfsk_setup(double sample_rate)
{
    unsigned cs = chip_samples(sample_rate);
    if (cs == _bfsk_cs && _bfsk_sin_lo != NULL)
        return;

    free(_bfsk_sin_lo);
    free(_bfsk_sin_hi);
    _bfsk_sin_lo = malloc(cs * sizeof(double));
    _bfsk_sin_hi = malloc(cs * sizeof(double));
    assert(_bfsk_sin_lo != NULL && _bfsk_sin_hi != NULL);
    for (unsigned s = 0; s < cs; s++) {
        double t = (double)s / sample_rate;
        _bfsk_sin_lo[s] = sin(2.0 * M_PI * FREQ_LO * t);
        _bfsk_sin_hi[s] = sin(2.0 * M_PI * FREQ_HI * t);
    }
    _bfsk_cs = cs;
}

void md_steg_teardown(void)
{
    free(_bfsk_sin_lo);
    free(_bfsk_sin_hi);
    _bfsk_sin_lo = NULL;
    _bfsk_sin_hi = NULL;
    _bfsk_cs     = 0;
}

/** Encode one bit by adding a precomputed BFSK tone burst at the given chip. */
static void encode_one_bit(double *output, unsigned signal_len,
                           unsigned chip_idx, unsigned cs,
                           const double *carrier)
{
    unsigned start = chip_idx * cs;
    for (unsigned s = 0; s < cs && (start + s) < signal_len; s++)
        output[start + s] += TONE_AMP * carrier[s];
}

static unsigned freq_encode(const double *host, double *output,
                            unsigned signal_len, double sample_rate,
                            const unsigned char *data, unsigned data_len,
                            unsigned flags)
{
    unsigned cs = chip_samples(sample_rate);
    unsigned capacity = freq_capacity_cs(signal_len, cs);
    if (data_len == 0 || capacity == 0)
        return 0;
    if (data_len > capacity)
        data_len = capacity;

    /* Copy host to output. */
    memcpy(output, host, signal_len * sizeof(double));

    /* Ensure cached sine carrier tables are current. */
    _bfsk_setup(sample_rate);

    /* Encode 32-bit length header (bit 31 carries payload type flag). */
    unsigned header = data_len | flags;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = (header >> i) & 1;
        encode_one_bit(output, signal_len, i, cs,
                       bit ? _bfsk_sin_hi : _bfsk_sin_lo);
    }

    /* Encode data bytes. */
    for (unsigned b = 0; b < data_len; b++) {
        unsigned char ch = data[b];
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned chip = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = (ch >> bit_idx) & 1;
            encode_one_bit(output, signal_len, chip, cs,
                           bit ? _bfsk_sin_hi : _bfsk_sin_lo);
        }
    }
    return data_len;
}

/** Decode one bit by correlating a chip against precomputed BFSK carriers.
 *  Optionally returns the winning correlation magnitude via @p corr_out. */
static unsigned decode_one_bit(const double *stego, unsigned signal_len,
                               unsigned chip_idx, unsigned cs,
                               const double *sin_lo, const double *sin_hi,
                               double *corr_out)
{
    unsigned start = chip_idx * cs;
    double corr_lo = 0.0, corr_hi = 0.0;
    for (unsigned s = 0; s < cs && (start + s) < signal_len; s++) {
        corr_lo += stego[start + s] * sin_lo[s];
        corr_hi += stego[start + s] * sin_hi[s];
    }
    if (fabs(corr_hi) > fabs(corr_lo)) {
        if (corr_out) *corr_out = fabs(corr_hi);
        return 1u;
    }
    if (corr_out) *corr_out = fabs(corr_lo);
    return 0u;
}

static unsigned freq_decode(const double *stego, unsigned signal_len,
                            double sample_rate,
                            unsigned char *data_out, unsigned max_len)
{
    unsigned cs = chip_samples(sample_rate);
    if (cs == 0 || max_len == 0)
        return 0;

    unsigned total_chips = signal_len / cs;
    if (total_chips <= HEADER_BITS)
        return 0;

    /* Ensure cached sine carrier tables are current. */
    _bfsk_setup(sample_rate);

    /* Read the 32-bit length header (bit 31 is payload type flag). */
    unsigned raw_header = 0;
    for (unsigned i = 0; i < HEADER_BITS; i++) {
        unsigned bit = decode_one_bit(stego, signal_len, i, cs,
                                      _bfsk_sin_lo, _bfsk_sin_hi,
                                      NULL);
        raw_header |= (bit << i);
    }
    unsigned msg_len = raw_header & LENGTH_MASK;

    /* Sanity check. */
    unsigned capacity = freq_capacity_cs(signal_len, cs);
    if (msg_len == 0 || msg_len > capacity)
        return 0;

    unsigned decode_len = msg_len;
    if (decode_len > max_len)
        decode_len = max_len;

    for (unsigned b = 0; b < decode_len; b++) {
        unsigned char ch = 0;
        for (unsigned bit_idx = 0; bit_idx < 8; bit_idx++) {
            unsigned chip = HEADER_BITS + b * 8 + bit_idx;
            unsigned bit = decode_one_bit(stego, signal_len, chip, cs,
                                          _bfsk_sin_lo, _bfsk_sin_hi,
                                          NULL);
            ch |= (unsigned char)(bit << bit_idx);
        }
        data_out[b] = ch;
    }

    return decode_len;
}

/* -----------------------------------------------------------------------
 * Spectrogram text steganography (hybrid LSB + visual)
 *
 * Combines LSB encoding (for machine-readable decode) with spectrogram
 * text art in the 18-24 kHz ultrasonic band (for visual verification).
 * Auto-upsamples to 48 kHz so Nyquist >= 24 kHz.
 *
 * Each bitmap column occupies a fixed 30 ms of audio; each character
 * is 8 columns (5 data + 3 spacing) = 240 ms.  The spectrogram text
 * signal is scaled to SPECTEXT_AMP (0.02) after generation.
 *
 * Bit layout: identical to LSB (32-bit header + payload).
 * ----------------------------------------------------------------------- */

#define SPECTEXT_FREQ_LO   18000.0   /**< Bottom of visual band (Hz). */
#define SPECTEXT_FREQ_HI   23500.0   /**< Top of visual band (Hz, below Nyquist). */
#define SPECTEXT_COL_MS    30.0      /**< Column duration (ms). */
#define SPECTEXT_COLS_PER_CHAR 8     /**< Bitmap columns per character. */
#define SPECTEXT_TARGET_SR 48000.0   /**< Output sample rate (Hz). */
#define SPECTEXT_AMP       0.02     /**< Spectrogram art amplitude. */
#define SPECTEXT_NORM_PEAK 0.9      /**< MD_spectrogram_text peak. */
#define SPECTEXT_PAD_SEC   0.25     /**< Silence before spectrogram text (s). */

/** Seconds per character in the spectrogram visual. */
#define SPECTEXT_SEC_PER_CHAR \
    (SPECTEXT_COL_MS / 1000.0 * SPECTEXT_COLS_PER_CHAR)

/** Maximum number of visually displayable characters for a given duration.
 *  Subtracts the leading pad from the available time. */
static unsigned spectext_vis_capacity(double duration_sec)
{
    double available = duration_sec - SPECTEXT_PAD_SEC;
    if (available <= 0.0) return 0;
    return (unsigned)(available / SPECTEXT_SEC_PER_CHAR);
}

/** Compute the output signal length at 48 kHz for a given input. */
static unsigned spectext_output_len(unsigned signal_len, double sample_rate)
{
    if (sample_rate >= SPECTEXT_TARGET_SR)
        return signal_len;
    return MD_resample_output_len(signal_len, sample_rate, SPECTEXT_TARGET_SR);
}

static unsigned spectext_capacity(unsigned signal_len, double sample_rate)
{
    double duration_sec = (double)signal_len / sample_rate;
    unsigned vis_cap = spectext_vis_capacity(duration_sec);
    unsigned out_len = spectext_output_len(signal_len, sample_rate);
    unsigned lsb_cap = lsb_capacity(out_len);
    return (vis_cap < lsb_cap) ? vis_cap : lsb_cap;
}

static unsigned spectext_encode(const double *host, double *output,
                                unsigned signal_len, double sample_rate,
                                const unsigned char *data, unsigned data_len,
                                unsigned flags)
{
    double duration_sec = (double)signal_len / sample_rate;
    unsigned out_len = spectext_output_len(signal_len, sample_rate);

    /* Step 1: Upsample host to 48 kHz if needed.
     * We need a mutable copy at 48 kHz to mix spectrogram art into
     * BEFORE LSB encoding (LSB must be the last step, because adding
     * any signal afterwards would disturb the LSB bits). */
    double *mixed = malloc(out_len * sizeof(double));
    assert(mixed != NULL);

    if (sample_rate < SPECTEXT_TARGET_SR) {
        MD_resample(host, signal_len, mixed, out_len,
                    sample_rate, SPECTEXT_TARGET_SR, 128, 10.0);
        /* Brickwall lowpass at the original Nyquist to eliminate any
         * residual spectral images from the resampler transition band. */
        MD_lowpass_brickwall(mixed, out_len,
                             sample_rate / 2.0, SPECTEXT_TARGET_SR);
    } else {
        memcpy(mixed, host, out_len * sizeof(double));
    }

    /* Step 2: Generate spectrogram text art and mix into host
     * (before LSB so the LSB bits remain undisturbed). */
    unsigned vis_chars = spectext_vis_capacity(duration_sec);
    if (vis_chars > 0) {
        /* For text payloads, show the message.
         * For binary payloads, show "[BIN <N>B]". */
        char vis_label[64];
        const char *vis_text;

        if (flags & PAYLOAD_TYPE_BIT) {
            snprintf(vis_label, sizeof(vis_label), "[BIN %uB]", data_len);
            vis_text = vis_label;
        } else {
            vis_text = (const char *)data;
        }

        unsigned text_len = (unsigned)strlen(vis_text);
        if (text_len > vis_chars)
            text_len = vis_chars;

        /* Build a null-terminated substring if truncated. */
        char *vis_substr = malloc(text_len + 1);
        assert(vis_substr != NULL);
        memcpy(vis_substr, vis_text, text_len);
        vis_substr[text_len] = '\0';

        double vis_duration = (double)text_len * SPECTEXT_SEC_PER_CHAR;

        double *specbuf = calloc(out_len, sizeof(double));
        assert(specbuf != NULL);

        MD_spectrogram_text(specbuf, out_len, vis_substr,
                            SPECTEXT_FREQ_LO, SPECTEXT_FREQ_HI,
                            vis_duration, SPECTEXT_TARGET_SR);

        /* Scale to steganographic amplitude and mix into host,
         * offset by the leading pad so text doesn't start at t=0. */
        double scale = SPECTEXT_AMP / SPECTEXT_NORM_PEAK;
        unsigned pad_samples = (unsigned)(SPECTEXT_PAD_SEC * SPECTEXT_TARGET_SR);
        unsigned spec_samples = (unsigned)(vis_duration * SPECTEXT_TARGET_SR);
        if (pad_samples + spec_samples > out_len)
            spec_samples = (out_len > pad_samples) ? out_len - pad_samples : 0;
        for (unsigned i = 0; i < spec_samples; i++)
            mixed[pad_samples + i] += specbuf[i] * scale;

        free(specbuf);
        free(vis_substr);
    }

    /* Step 3: LSB-encode the full message (last step — must not
     * add any signal afterwards or the LSB bits get disturbed). */
    unsigned encoded = lsb_encode(mixed, output, out_len,
                                   data, data_len, flags);
    free(mixed);

    return encoded;
}

/** Probe for ultrasonic energy in the 18-24 kHz band.
 *  Computes RMS of the signal bandpass-filtered to that range using
 *  a simple DFT bin summation approach. */
static double spectext_ultrasonic_rms(const double *signal,
                                       unsigned signal_len,
                                       double sample_rate)
{
    /* Use a short window from the signal to check for ultrasonic energy.
     * We compute a DFT and sum energy in the 18-24 kHz bins. */
    unsigned fft_len = 4096;
    if (fft_len > signal_len)
        fft_len = signal_len;

    /* Compute magnitude spectrum (3-arg: signal, N, mag_out). */
    unsigned spec_len = fft_len / 2 + 1;
    double *mag = malloc(spec_len * sizeof(double));
    assert(mag != NULL);

    /* Use a segment from the middle of the signal. */
    unsigned offset = 0;
    if (signal_len > fft_len)
        offset = (signal_len - fft_len) / 2;

    double *windowed = malloc(fft_len * sizeof(double));
    double *win = malloc(fft_len * sizeof(double));
    assert(windowed != NULL && win != NULL);
    MD_Gen_Hann_Win(win, fft_len);
    for (unsigned i = 0; i < fft_len; i++)
        windowed[i] = signal[offset + i] * win[i];
    free(win);

    MD_magnitude_spectrum(windowed, fft_len, mag);

    /* Sum energy in the 18-24 kHz range. */
    double bin_hz = sample_rate / (double)fft_len;
    unsigned bin_lo = (unsigned)(SPECTEXT_FREQ_LO / bin_hz);
    unsigned bin_hi = (unsigned)(SPECTEXT_FREQ_HI / bin_hz);
    if (bin_hi >= spec_len)
        bin_hi = spec_len - 1;

    double energy = 0.0;
    unsigned count = 0;
    for (unsigned k = bin_lo; k <= bin_hi; k++) {
        energy += mag[k] * mag[k];
        count++;
    }

    free(windowed);
    free(mag);

    if (count == 0) return 0.0;
    return sqrt(energy / (double)count);
}

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

unsigned MD_steg_capacity(unsigned signal_len, double sample_rate, int method)
{
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND ||
           method == MD_STEG_SPECTEXT);

    if (method == MD_STEG_LSB)
        return lsb_capacity(signal_len);
    else if (method == MD_STEG_SPECTEXT)
        return spectext_capacity(signal_len, sample_rate);
    else
        return freq_capacity(signal_len, sample_rate);
}

/** Shared encode logic for both text and binary public API functions. */
static unsigned encode_common(const double *host, double *output,
                              unsigned signal_len, double sample_rate,
                              const unsigned char *data, unsigned data_len,
                              int method, unsigned flags)
{
    assert(host != NULL);
    assert(output != NULL);
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(data != NULL);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND ||
           method == MD_STEG_SPECTEXT);

    if (method == MD_STEG_FREQ_BAND)
        assert(sample_rate >= 40000.0 &&
               "frequency-band steganography requires sample_rate >= 40 kHz");

    if (method == MD_STEG_LSB)
        return lsb_encode(host, output, signal_len, data, data_len, flags);
    else if (method == MD_STEG_SPECTEXT)
        return spectext_encode(host, output, signal_len, sample_rate,
                               data, data_len, flags);
    else
        return freq_encode(host, output, signal_len, sample_rate,
                           data, data_len, flags);
}

unsigned MD_steg_encode_bytes(const double *host, double *output,
                              unsigned signal_len, double sample_rate,
                              const unsigned char *data, unsigned data_len,
                              int method)
{
    return encode_common(host, output, signal_len, sample_rate,
                         data, data_len, method, PAYLOAD_TYPE_BIT);
}

unsigned MD_steg_decode_bytes(const double *stego, unsigned signal_len,
                              double sample_rate,
                              unsigned char *data_out, unsigned max_len,
                              int method)
{
    assert(stego != NULL);
    assert(signal_len > 0);
    assert(sample_rate > 0.0);
    assert(data_out != NULL);
    assert(max_len > 0);
    assert(method == MD_STEG_LSB || method == MD_STEG_FREQ_BAND ||
           method == MD_STEG_SPECTEXT);

    if (method == MD_STEG_LSB || method == MD_STEG_SPECTEXT)
        return lsb_decode(stego, signal_len, data_out, max_len);
    else
        return freq_decode(stego, signal_len, sample_rate,
                           data_out, max_len);
}

unsigned MD_steg_encode(const double *host, double *output,
                        unsigned signal_len, double sample_rate,
                        const char *message, int method)
{
    assert(message != NULL);
    return encode_common(host, output, signal_len, sample_rate,
                         (const unsigned char *)message,
                         (unsigned)strlen(message), method, 0);
}

unsigned MD_steg_decode(const double *stego, unsigned signal_len,
                        double sample_rate,
                        char *message_out, unsigned max_msg_len,
                        int method)
{
    assert(message_out != NULL);
    assert(max_msg_len > 0);
    unsigned decoded = MD_steg_decode_bytes(stego, signal_len, sample_rate,
                                            (unsigned char *)message_out,
                                            max_msg_len - 1, method);
    message_out[decoded] = '\0';
    return decoded;
}

int MD_steg_detect(const double *signal, unsigned signal_len,
                   double sample_rate, int *payload_type_out)
{
    assert(signal != NULL);
    assert(signal_len > 0);
    assert(sample_rate > 0.0);

    int found_method = -1;
    unsigned found_header = 0;

    /* --- BFSK probe (only when sample_rate >= 40 kHz) --- */
    if (sample_rate >= 40000.0) {
        unsigned cs = chip_samples(sample_rate);
        if (cs > 0) {
            unsigned total_chips = signal_len / cs;
            if (total_chips > HEADER_BITS) {
                _bfsk_setup(sample_rate);

                unsigned raw_header = 0;
                double corr_sum = 0.0;
                for (unsigned i = 0; i < HEADER_BITS; i++) {
                    double corr;
                    unsigned bit = decode_one_bit(
                        signal, signal_len, i, cs,
                        _bfsk_sin_lo, _bfsk_sin_hi, &corr);
                    raw_header |= (bit << i);
                    corr_sum += corr;
                }
                unsigned msg_len = raw_header & LENGTH_MASK;
                unsigned capacity = freq_capacity_cs(signal_len, cs);
                double avg_corr = corr_sum / HEADER_BITS;
                double threshold = 0.25 * TONE_AMP * (double)cs / 2.0;

                if (msg_len > 0 && msg_len <= capacity &&
                    avg_corr >= threshold) {
                    found_method = MD_STEG_FREQ_BAND;
                    found_header = raw_header;
                }
            }
        }
    }

    /* --- LSB probe (also covers spectext, which uses LSB + ultrasonic art) --- */
    if (found_method < 0 && signal_len > HEADER_BITS) {
        unsigned raw_header = lsb_read_header(signal);
        unsigned msg_len = raw_header & LENGTH_MASK;
        unsigned capacity = lsb_capacity(signal_len);

        if (msg_len > 0 && msg_len <= capacity) {
            /* Check for ultrasonic energy to distinguish spectext from LSB. */
            if (sample_rate >= SPECTEXT_TARGET_SR) {
                double ultra_rms = spectext_ultrasonic_rms(signal, signal_len,
                                                           sample_rate);
                if (ultra_rms > 1e-4) {
                    found_method = MD_STEG_SPECTEXT;
                    found_header = raw_header;
                } else {
                    found_method = MD_STEG_LSB;
                    found_header = raw_header;
                }
            } else {
                found_method = MD_STEG_LSB;
                found_header = raw_header;
            }
        }
    }

    if (found_method >= 0 && payload_type_out != NULL)
        *payload_type_out = (found_header >> 31) ? MD_STEG_TYPE_BINARY
                                                 : MD_STEG_TYPE_TEXT;

    return found_method;
}
