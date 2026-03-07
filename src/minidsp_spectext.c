/**
 * @file minidsp_spectext.c
 * @brief Spectrogram text art: synthesise audio that displays readable text
 *        when viewed as a spectrogram.
 *
 * A 5x7 bitmap font rasterises each printable ASCII character.  Each bitmap
 * column becomes a time slice; each "on" pixel becomes a sine wave at the
 * corresponding frequency.  A raised-cosine crossfade suppresses clicks at
 * column boundaries.
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * 5x7 Bitmap Font (printable ASCII 32–126)
 *
 * Each character occupies 5 columns.  Each column is one byte: bit 0 is
 * the top row, bit 6 is the bottom row.  Bits 7 is unused.
 * -----------------------------------------------------------------------*/

static const unsigned char font_5x7[95][5] = {
    /* 32 ' ' */ {0x00, 0x00, 0x00, 0x00, 0x00},
    /* 33 '!' */ {0x00, 0x00, 0x5F, 0x00, 0x00},
    /* 34 '"' */ {0x00, 0x07, 0x00, 0x07, 0x00},
    /* 35 '#' */ {0x14, 0x7F, 0x14, 0x7F, 0x14},
    /* 36 '$' */ {0x24, 0x2A, 0x7F, 0x2A, 0x12},
    /* 37 '%' */ {0x23, 0x13, 0x08, 0x64, 0x62},
    /* 38 '&' */ {0x36, 0x49, 0x55, 0x22, 0x50},
    /* 39 ''' */ {0x00, 0x05, 0x03, 0x00, 0x00},
    /* 40 '(' */ {0x00, 0x1C, 0x22, 0x41, 0x00},
    /* 41 ')' */ {0x00, 0x41, 0x22, 0x1C, 0x00},
    /* 42 '*' */ {0x14, 0x08, 0x3E, 0x08, 0x14},
    /* 43 '+' */ {0x08, 0x08, 0x3E, 0x08, 0x08},
    /* 44 ',' */ {0x00, 0x50, 0x30, 0x00, 0x00},
    /* 45 '-' */ {0x08, 0x08, 0x08, 0x08, 0x08},
    /* 46 '.' */ {0x00, 0x60, 0x60, 0x00, 0x00},
    /* 47 '/' */ {0x20, 0x10, 0x08, 0x04, 0x02},
    /* 48 '0' */ {0x3E, 0x51, 0x49, 0x45, 0x3E},
    /* 49 '1' */ {0x00, 0x42, 0x7F, 0x40, 0x00},
    /* 50 '2' */ {0x42, 0x61, 0x51, 0x49, 0x46},
    /* 51 '3' */ {0x21, 0x41, 0x45, 0x4B, 0x31},
    /* 52 '4' */ {0x18, 0x14, 0x12, 0x7F, 0x10},
    /* 53 '5' */ {0x27, 0x45, 0x45, 0x45, 0x39},
    /* 54 '6' */ {0x3C, 0x4A, 0x49, 0x49, 0x30},
    /* 55 '7' */ {0x01, 0x71, 0x09, 0x05, 0x03},
    /* 56 '8' */ {0x36, 0x49, 0x49, 0x49, 0x36},
    /* 57 '9' */ {0x06, 0x49, 0x49, 0x29, 0x1E},
    /* 58 ':' */ {0x00, 0x36, 0x36, 0x00, 0x00},
    /* 59 ';' */ {0x00, 0x56, 0x36, 0x00, 0x00},
    /* 60 '<' */ {0x08, 0x14, 0x22, 0x41, 0x00},
    /* 61 '=' */ {0x14, 0x14, 0x14, 0x14, 0x14},
    /* 62 '>' */ {0x00, 0x41, 0x22, 0x14, 0x08},
    /* 63 '?' */ {0x02, 0x01, 0x51, 0x09, 0x06},
    /* 64 '@' */ {0x32, 0x49, 0x79, 0x41, 0x3E},
    /* 65 'A' */ {0x7E, 0x11, 0x11, 0x11, 0x7E},
    /* 66 'B' */ {0x7F, 0x49, 0x49, 0x49, 0x36},
    /* 67 'C' */ {0x3E, 0x41, 0x41, 0x41, 0x22},
    /* 68 'D' */ {0x7F, 0x41, 0x41, 0x22, 0x1C},
    /* 69 'E' */ {0x7F, 0x49, 0x49, 0x49, 0x41},
    /* 70 'F' */ {0x7F, 0x09, 0x09, 0x09, 0x01},
    /* 71 'G' */ {0x3E, 0x41, 0x49, 0x49, 0x7A},
    /* 72 'H' */ {0x7F, 0x08, 0x08, 0x08, 0x7F},
    /* 73 'I' */ {0x00, 0x41, 0x7F, 0x41, 0x00},
    /* 74 'J' */ {0x20, 0x40, 0x41, 0x3F, 0x01},
    /* 75 'K' */ {0x7F, 0x08, 0x14, 0x22, 0x41},
    /* 76 'L' */ {0x7F, 0x40, 0x40, 0x40, 0x40},
    /* 77 'M' */ {0x7F, 0x02, 0x0C, 0x02, 0x7F},
    /* 78 'N' */ {0x7F, 0x04, 0x08, 0x10, 0x7F},
    /* 79 'O' */ {0x3E, 0x41, 0x41, 0x41, 0x3E},
    /* 80 'P' */ {0x7F, 0x09, 0x09, 0x09, 0x06},
    /* 81 'Q' */ {0x3E, 0x41, 0x51, 0x21, 0x5E},
    /* 82 'R' */ {0x7F, 0x09, 0x19, 0x29, 0x46},
    /* 83 'S' */ {0x46, 0x49, 0x49, 0x49, 0x31},
    /* 84 'T' */ {0x01, 0x01, 0x7F, 0x01, 0x01},
    /* 85 'U' */ {0x3F, 0x40, 0x40, 0x40, 0x3F},
    /* 86 'V' */ {0x1F, 0x20, 0x40, 0x20, 0x1F},
    /* 87 'W' */ {0x3F, 0x40, 0x38, 0x40, 0x3F},
    /* 88 'X' */ {0x63, 0x14, 0x08, 0x14, 0x63},
    /* 89 'Y' */ {0x07, 0x08, 0x70, 0x08, 0x07},
    /* 90 'Z' */ {0x61, 0x51, 0x49, 0x45, 0x43},
    /* 91 '[' */ {0x00, 0x7F, 0x41, 0x41, 0x00},
    /* 92 '\' */ {0x02, 0x04, 0x08, 0x10, 0x20},
    /* 93 ']' */ {0x00, 0x41, 0x41, 0x7F, 0x00},
    /* 94 '^' */ {0x04, 0x02, 0x01, 0x02, 0x04},
    /* 95 '_' */ {0x40, 0x40, 0x40, 0x40, 0x40},
    /* 96 '`' */ {0x00, 0x01, 0x02, 0x04, 0x00},
    /* 97 'a' */ {0x20, 0x54, 0x54, 0x54, 0x78},
    /* 98 'b' */ {0x7F, 0x48, 0x44, 0x44, 0x38},
    /* 99 'c' */ {0x38, 0x44, 0x44, 0x44, 0x20},
    /*100 'd' */ {0x38, 0x44, 0x44, 0x48, 0x7F},
    /*101 'e' */ {0x38, 0x54, 0x54, 0x54, 0x18},
    /*102 'f' */ {0x08, 0x7E, 0x09, 0x01, 0x02},
    /*103 'g' */ {0x0C, 0x52, 0x52, 0x52, 0x3E},
    /*104 'h' */ {0x7F, 0x08, 0x04, 0x04, 0x78},
    /*105 'i' */ {0x00, 0x44, 0x7D, 0x40, 0x00},
    /*106 'j' */ {0x20, 0x40, 0x44, 0x3D, 0x00},
    /*107 'k' */ {0x7F, 0x10, 0x28, 0x44, 0x00},
    /*108 'l' */ {0x00, 0x41, 0x7F, 0x40, 0x00},
    /*109 'm' */ {0x7C, 0x04, 0x18, 0x04, 0x78},
    /*110 'n' */ {0x7C, 0x08, 0x04, 0x04, 0x78},
    /*111 'o' */ {0x38, 0x44, 0x44, 0x44, 0x38},
    /*112 'p' */ {0x7C, 0x14, 0x14, 0x14, 0x08},
    /*113 'q' */ {0x08, 0x14, 0x14, 0x18, 0x7C},
    /*114 'r' */ {0x7C, 0x08, 0x04, 0x04, 0x08},
    /*115 's' */ {0x48, 0x54, 0x54, 0x54, 0x20},
    /*116 't' */ {0x04, 0x3F, 0x44, 0x40, 0x20},
    /*117 'u' */ {0x3C, 0x40, 0x40, 0x20, 0x7C},
    /*118 'v' */ {0x1C, 0x20, 0x40, 0x20, 0x1C},
    /*119 'w' */ {0x3C, 0x40, 0x30, 0x40, 0x3C},
    /*120 'x' */ {0x44, 0x28, 0x10, 0x28, 0x44},
    /*121 'y' */ {0x0C, 0x50, 0x50, 0x50, 0x3C},
    /*122 'z' */ {0x44, 0x64, 0x54, 0x4C, 0x44},
    /*123 '{' */ {0x00, 0x08, 0x36, 0x41, 0x00},
    /*124 '|' */ {0x00, 0x00, 0x7F, 0x00, 0x00},
    /*125 '}' */ {0x00, 0x41, 0x36, 0x08, 0x00},
    /*126 '~' */ {0x10, 0x08, 0x08, 0x10, 0x10},
};

/** Return the pixel state (0 or 1) at bitmap coordinate (col, row)
 *  for the given string.  Characters are 5 columns wide with 3-column
 *  spacing.  Rows are numbered 0 (top) to 6 (bottom). */
static int pixel_at(const char *text, unsigned len, unsigned col, unsigned row)
{
    unsigned char_width = 8;  /* 5 data columns + 3 spacing columns */
    unsigned char_idx = col / char_width;
    unsigned col_in_char = col % char_width;

    if (char_idx >= len) return 0;
    if (col_in_char >= 5) return 0;  /* spacing column */

    char ch = text[char_idx];
    if (ch < 32 || ch > 126) ch = '?';  /* replace unprintable */

    unsigned glyph_idx = (unsigned)(ch - 32);
    unsigned char col_byte = font_5x7[glyph_idx][col_in_char];
    return (col_byte >> row) & 1;
}

/* -----------------------------------------------------------------------
 * Public API: spectrogram text synthesis
 * -----------------------------------------------------------------------*/

unsigned MD_spectrogram_text(double *output, unsigned max_len,
                             const char *text,
                             double freq_lo, double freq_hi,
                             double duration_sec, double sample_rate)
{
    assert(output != NULL);
    assert(text != NULL);
    assert(text[0] != '\0');
    assert(freq_lo > 0.0);
    assert(freq_lo < freq_hi);
    assert(freq_hi < sample_rate / 2.0);
    assert(duration_sec > 0.0);
    assert(sample_rate > 0.0);

    unsigned len = (unsigned)strlen(text);
    unsigned grid_cols = len * 8 - 3;  /* 5 data + 3 space per char, minus trailing spaces */
    unsigned grid_rows = 7;

    unsigned col_samples = (unsigned)(duration_sec / (double)grid_cols * sample_rate);
    if (col_samples < 1) col_samples = 1;

    unsigned total_samples = col_samples * grid_cols;
    assert(max_len >= total_samples);

    /* Crossfade: 3 ms raised-cosine, clamped to half a column */
    unsigned fade_samples = (unsigned)(0.003 * sample_rate);
    if (fade_samples > col_samples / 2) fade_samples = col_samples / 2;

    /* Per-row frequency: row 0 → freq_hi, row 6 → freq_lo */
    double row_freq[7];
    for (unsigned r = 0; r < grid_rows; r++) {
        row_freq[r] = freq_hi - (double)r / (double)(grid_rows - 1)
                       * (freq_hi - freq_lo);
    }

    /* Per-row phase accumulator (for continuous tones) */
    double phase[7] = {0};

    /* Zero the output buffer */
    memset(output, 0, total_samples * sizeof(double));

    for (unsigned c = 0; c < grid_cols; c++) {
        unsigned base = c * col_samples;

        for (unsigned r = 0; r < grid_rows; r++) {
            int on_now  = pixel_at(text, len, c, r);
            int on_prev = (c > 0) ? pixel_at(text, len, c - 1, r) : 0;
            int on_next = (c + 1 < grid_cols) ? pixel_at(text, len, c + 1, r) : 0;

            double dp = 2.0 * M_PI * row_freq[r] / sample_rate;

            for (unsigned s = 0; s < col_samples; s++) {
                if (!on_now) {
                    /* Advance phase even when off so re-entry is smooth */
                    phase[r] = fmod(phase[r] + dp, 2.0 * M_PI);
                    continue;
                }

                double envelope = 1.0;

                /* Fade-in at column start if previous column was off */
                if (!on_prev && s < fade_samples && fade_samples > 0) {
                    double t = (double)s / (double)fade_samples;
                    envelope *= 0.5 * (1.0 - cos(M_PI * t));
                }

                /* Fade-out at column end if next column is off */
                if (!on_next && s >= col_samples - fade_samples && fade_samples > 0) {
                    double t = (double)(col_samples - 1 - s) / (double)fade_samples;
                    envelope *= 0.5 * (1.0 - cos(M_PI * t));
                }

                output[base + s] += envelope * sin(phase[r]);
                phase[r] = fmod(phase[r] + dp, 2.0 * M_PI);
            }
        }
    }

    /* Normalize to 0.9 peak amplitude */
    double peak = 0.0;
    for (unsigned i = 0; i < total_samples; i++) {
        double a = fabs(output[i]);
        if (a > peak) peak = a;
    }
    if (peak > 0.0) {
        double scale = 0.9 / peak;
        for (unsigned i = 0; i < total_samples; i++) {
            output[i] *= scale;
        }
    }

    return total_samples;
}
