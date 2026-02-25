/**
 * @file minidsp_dtmf.c
 * @brief DTMF tone detection and generation.
 *
 * Detection uses frame-based FFT magnitude analysis with a state machine
 * enforcing ITU-T Q.24 minimum timing (40 ms tone-on, 40 ms pause).
 *
 * Generation superimposes row + column sinusoids via MD_sine_wave().
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * DTMF frequency tables
 * ----------------------------------------------------------------------- */

static const double dtmf_row_freqs[4] = {697.0, 770.0, 852.0, 941.0};
static const double dtmf_col_freqs[4] = {1209.0, 1336.0, 1477.0, 1633.0};

/*  Keypad layout (row, col):
 *      1209   1336   1477   1633
 * 697:   1      2      3      A
 * 770:   4      5      6      B
 * 852:   7      8      9      C
 * 941:   *      0      #      D
 */
static const char dtmf_keypad[4][4] = {
    {'1', '2', '3', 'A'},
    {'4', '5', '6', 'B'},
    {'7', '8', '9', 'C'},
    {'*', '0', '#', 'D'}
};

/* -----------------------------------------------------------------------
 * Internal helpers
 * ----------------------------------------------------------------------- */

/** Map a DTMF character to its row and column frequencies.
 *  Asserts on invalid characters (structural misuse). */
static void dtmf_char_to_freqs(char ch, double *row_freq, double *col_freq)
{
    int row = -1, col = -1;
    switch (ch) {
    case '1': row = 0; col = 0; break;
    case '2': row = 0; col = 1; break;
    case '3': row = 0; col = 2; break;
    case 'A': case 'a': row = 0; col = 3; break;
    case '4': row = 1; col = 0; break;
    case '5': row = 1; col = 1; break;
    case '6': row = 1; col = 2; break;
    case 'B': case 'b': row = 1; col = 3; break;
    case '7': row = 2; col = 0; break;
    case '8': row = 2; col = 1; break;
    case '9': row = 2; col = 2; break;
    case 'C': case 'c': row = 2; col = 3; break;
    case '*': row = 3; col = 0; break;
    case '0': row = 3; col = 1; break;
    case '#': row = 3; col = 2; break;
    case 'D': case 'd': row = 3; col = 3; break;
    default:
        assert(0 && "invalid DTMF character");
    }
    *row_freq = dtmf_row_freqs[row];
    *col_freq = dtmf_col_freqs[col];
}

/** Peak magnitude in bins [bin-1, bin, bin+1], clamped to [0, num_bins). */
static double peak_near_bin(const double *mag, unsigned num_bins, unsigned bin)
{
    double peak = mag[bin];
    if (bin > 0 && mag[bin - 1] > peak)
        peak = mag[bin - 1];
    if (bin + 1 < num_bins && mag[bin + 1] > peak)
        peak = mag[bin + 1];
    return peak;
}

/** Detect the DTMF digit present in a single normalised magnitude frame.
 *  Returns the digit character, or '\0' if no valid DTMF pair is found. */
static char detect_frame(const double *mag, unsigned num_bins,
                         unsigned N, double sample_rate)
{
    /* Compute mean magnitude (exclude DC and Nyquist) for threshold. */
    double sum = 0.0;
    for (unsigned k = 1; k + 1 < num_bins; k++)
        sum += mag[k];
    double mean_mag = (num_bins > 2) ? sum / (double)(num_bins - 2) : 0.0;
    double threshold = mean_mag * 8.0;   /* ~18 dB above mean */

    /* Measure magnitude at each DTMF row frequency. */
    double row_mags[4];
    for (int r = 0; r < 4; r++) {
        unsigned bin = (unsigned)(dtmf_row_freqs[r] * N / sample_rate + 0.5);
        if (bin >= num_bins) bin = num_bins - 1;
        row_mags[r] = peak_near_bin(mag, num_bins, bin);
    }

    /* Measure magnitude at each DTMF column frequency. */
    double col_mags[4];
    for (int c = 0; c < 4; c++) {
        unsigned bin = (unsigned)(dtmf_col_freqs[c] * N / sample_rate + 0.5);
        if (bin >= num_bins) bin = num_bins - 1;
        col_mags[c] = peak_near_bin(mag, num_bins, bin);
    }

    /* Find the strongest row and column that exceed the threshold. */
    int best_row = -1;
    double best_row_mag = 0.0;
    for (int r = 0; r < 4; r++) {
        if (row_mags[r] > threshold && row_mags[r] > best_row_mag) {
            best_row = r;
            best_row_mag = row_mags[r];
        }
    }

    int best_col = -1;
    double best_col_mag = 0.0;
    for (int c = 0; c < 4; c++) {
        if (col_mags[c] > threshold && col_mags[c] > best_col_mag) {
            best_col = c;
            best_col_mag = col_mags[c];
        }
    }

    if (best_row < 0 || best_col < 0)
        return '\0';

    return dtmf_keypad[best_row][best_col];
}

/* -----------------------------------------------------------------------
 * Public API
 * ----------------------------------------------------------------------- */

unsigned MD_dtmf_detect(const double *signal, unsigned signal_len,
                        double sample_rate,
                        MD_DTMFTone *tones_out, unsigned max_tones)
{
    assert(signal);
    assert(tones_out);
    assert(signal_len > 0);
    assert(sample_rate >= 4000.0);
    assert(max_tones > 0);

    /* Pick FFT size: need enough resolution to separate DTMF row
     * pairs (73 Hz minimum gap → need < 37 Hz resolution) but the
     * window must stay shorter than the 40 ms Q.24 minimum pause so
     * that the frame-based state machine can resolve inter-digit gaps.
     * Target: largest power of two with window <= 35 ms. */
    unsigned max_n = (unsigned)(0.035 * sample_rate);
    unsigned N = 128;
    while (N * 2 <= max_n) N <<= 1;
    unsigned hop = N / 4;
    unsigned num_bins = N / 2 + 1;

    unsigned num_frames = (signal_len >= N)
                        ? (signal_len - N) / hop + 1
                        : 0;
    if (num_frames == 0)
        return 0;

    /* ITU-T Q.24: 40 ms minimum tone-on and inter-digit pause.
     * Compute the number of frames whose analysis window fits entirely
     * within a 40 ms interval, with a floor of 2 to debounce noise. */
    unsigned q24_samples = (unsigned)(0.040 * sample_rate);
    unsigned min_on_frames  = (q24_samples >= N)
                            ? (q24_samples - N) / hop + 1 : 1;
    if (min_on_frames < 2) min_on_frames = 2;
    unsigned min_off_frames = min_on_frames;

    /* Working buffers. */
    double *window = malloc(N * sizeof(double));
    double *frame  = malloc(N * sizeof(double));
    double *mag    = malloc(num_bins * sizeof(double));
    assert(window && frame && mag);

    MD_Gen_Hann_Win(window, N);

    /* State machine. */
    enum { IDLE, PENDING, ACTIVE } state = IDLE;
    char   current_digit    = '\0';
    unsigned on_count        = 0;
    unsigned off_count       = 0;
    unsigned tone_start_frame = 0;
    unsigned tone_end_frame   = 0;
    unsigned num_tones        = 0;

    for (unsigned f = 0; f < num_frames && num_tones < max_tones; f++) {
        unsigned start = f * hop;

        /* Window the frame. */
        for (unsigned i = 0; i < N; i++)
            frame[i] = signal[start + i] * window[i];

        /* Magnitude spectrum. */
        MD_magnitude_spectrum(frame, N, mag);

        /* Normalise to single-sided amplitude. */
        for (unsigned k = 0; k < num_bins; k++) {
            mag[k] /= (double)N;
            if (k > 0 && k < N / 2)
                mag[k] *= 2.0;
        }

        char digit = detect_frame(mag, num_bins, N, sample_rate);

        switch (state) {
        case IDLE:
            if (digit != '\0') {
                current_digit    = digit;
                on_count         = 1;
                tone_start_frame = f;
                state            = PENDING;
            }
            break;

        case PENDING:
            if (digit == current_digit) {
                on_count++;
                if (on_count >= min_on_frames) {
                    state          = ACTIVE;
                    tone_end_frame = f;
                    off_count      = 0;
                }
            } else if (digit != '\0') {
                /* Different digit — restart. */
                current_digit    = digit;
                on_count         = 1;
                tone_start_frame = f;
            } else {
                state         = IDLE;
                current_digit = '\0';
            }
            break;

        case ACTIVE:
            if (digit == current_digit && off_count == 0) {
                /* Same digit, no gap — tone continues. */
                tone_end_frame = f;
            } else if (digit == current_digit
                       && off_count >= min_off_frames) {
                /* Same digit reappeared after a gap >= Q.24 pause.
                 * This is a new instance of the same digit.
                 * Emit the current tone and start fresh. */
                tones_out[num_tones].digit   = current_digit;
                tones_out[num_tones].start_s =
                    (double)(tone_start_frame * hop) / sample_rate;
                tones_out[num_tones].end_s =
                    (double)((tone_end_frame + 1) * hop) / sample_rate;
                num_tones++;

                current_digit    = digit;
                on_count         = 1;
                tone_start_frame = f;
                state            = PENDING;
            } else if (digit == current_digit) {
                /* Brief interruption (< Q.24 pause), tolerate. */
                off_count      = 0;
                tone_end_frame = f;
            } else {
                off_count++;
                if (off_count >= min_off_frames) {
                    /* Emit the completed tone. */
                    tones_out[num_tones].digit   = current_digit;
                    tones_out[num_tones].start_s =
                        (double)(tone_start_frame * hop) / sample_rate;
                    tones_out[num_tones].end_s =
                        (double)((tone_end_frame + 1) * hop) / sample_rate;
                    num_tones++;

                    state         = IDLE;
                    current_digit = '\0';
                }
            }
            break;
        }
    }

    /* Emit a tone still active at end-of-signal.
     * PENDING is only emitted if it already meets the Q.24 minimum-on
     * requirement (on_count >= min_on_frames) to avoid false trailing
     * detections. */
    if (num_tones < max_tones) {
        if (state == ACTIVE) {
            tones_out[num_tones].digit   = current_digit;
            tones_out[num_tones].start_s =
                (double)(tone_start_frame * hop) / sample_rate;
            tones_out[num_tones].end_s =
                (double)((tone_end_frame + 1) * hop) / sample_rate;
            num_tones++;
        } else if (state == PENDING && on_count >= min_on_frames) {
            unsigned last = tone_start_frame + on_count - 1;
            tones_out[num_tones].digit   = current_digit;
            tones_out[num_tones].start_s =
                (double)(tone_start_frame * hop) / sample_rate;
            tones_out[num_tones].end_s =
                (double)((last + 1) * hop) / sample_rate;
            num_tones++;
        }
    }

    free(mag);
    free(frame);
    free(window);
    return num_tones;
}

void MD_dtmf_generate(double *output, const char *digits,
                      double sample_rate,
                      unsigned tone_ms, unsigned pause_ms)
{
    assert(output);
    assert(digits);
    assert(sample_rate > 0);
    assert(tone_ms >= 40);
    assert(pause_ms >= 40);

    unsigned num_digits   = (unsigned)strlen(digits);
    unsigned tone_samples = (unsigned)(tone_ms * sample_rate / 1000.0);
    unsigned pause_samples = (unsigned)(pause_ms * sample_rate / 1000.0);
    unsigned total = MD_dtmf_signal_length(num_digits, sample_rate,
                                           tone_ms, pause_ms);

    /* Zero-fill the entire output (silences between tones). */
    memset(output, 0, total * sizeof(double));

    /* Temporary buffer for the column sinusoid. */
    double *col_tone = malloc(tone_samples * sizeof(double));
    assert(col_tone);

    unsigned offset = 0;
    for (unsigned d = 0; d < num_digits; d++) {
        double row_freq, col_freq;
        dtmf_char_to_freqs(digits[d], &row_freq, &col_freq);

        /* Row tone directly into output at amplitude 0.5. */
        MD_sine_wave(output + offset, tone_samples, 0.5, row_freq, sample_rate);

        /* Column tone into temp, then add. */
        MD_sine_wave(col_tone, tone_samples, 0.5, col_freq, sample_rate);
        for (unsigned i = 0; i < tone_samples; i++)
            output[offset + i] += col_tone[i];

        offset += tone_samples + pause_samples;
    }

    free(col_tone);
}

unsigned MD_dtmf_signal_length(unsigned num_digits, double sample_rate,
                               unsigned tone_ms, unsigned pause_ms)
{
    assert(sample_rate > 0);
    if (num_digits == 0) return 0;

    unsigned tone_samples  = (unsigned)(tone_ms * sample_rate / 1000.0);
    unsigned pause_samples = (unsigned)(pause_ms * sample_rate / 1000.0);
    return num_digits * tone_samples
         + (num_digits - 1) * pause_samples;
}
