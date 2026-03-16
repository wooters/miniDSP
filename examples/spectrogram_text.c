/**
 * @file spectrogram_text.c
 * @brief Example: generate audio whose spectrogram displays readable text.
 *
 * Usage:
 *   ./spectrogram_text [TEXT] [--colormap hot|grayscale]
 *
 * Outputs:
 *   spectrogram_text.wav  — synthesised audio (WAV, 16-bit PCM)
 *   spectrogram_text.html — interactive Plotly spectrogram
 *   spectrogram_text.png  — static PNG spectrogram image
 *
 * Build and run (from the repository root):
 *   make -C examples spectrogram_text
 *   cd examples && ./spectrogram_text "HELLO"
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "fileio.h"
#include "plot_html.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../third_party/stb_image_write.h"

/* -----------------------------------------------------------------------
 * Hot 256-entry colormap (R, G, B): black → red → yellow → white
 * -----------------------------------------------------------------------*/

static const unsigned char hot_lut[256][3] = {
    {0,0,0},{3,0,0},{5,0,0},{8,0,0},{11,0,0},{13,0,0},{16,0,0},{19,0,0},
    {21,0,0},{24,0,0},{27,0,0},{30,0,0},{32,0,0},{35,0,0},{38,0,0},{40,0,0},
    {43,0,0},{46,0,0},{48,0,0},{51,0,0},{54,0,0},{56,0,0},{59,0,0},{62,0,0},
    {64,0,0},{67,0,0},{70,0,0},{72,0,0},{75,0,0},{78,0,0},{81,0,0},{83,0,0},
    {86,0,0},{89,0,0},{91,0,0},{94,0,0},{97,0,0},{99,0,0},{102,0,0},
    {105,0,0},{107,0,0},{110,0,0},{113,0,0},{115,0,0},{118,0,0},{121,0,0},
    {123,0,0},{126,0,0},{129,0,0},{132,0,0},{134,0,0},{137,0,0},{140,0,0},
    {142,0,0},{145,0,0},{148,0,0},{150,0,0},{153,0,0},{156,0,0},{158,0,0},
    {161,0,0},{164,0,0},{166,0,0},{169,0,0},{172,0,0},{174,0,0},{177,0,0},
    {180,0,0},{183,0,0},{185,0,0},{188,0,0},{191,0,0},{193,0,0},{196,0,0},
    {199,0,0},{201,0,0},{204,0,0},{207,0,0},{209,0,0},{212,0,0},{215,0,0},
    {217,0,0},{220,0,0},{223,0,0},{225,0,0},{228,0,0},{231,0,0},{234,0,0},
    {236,0,0},{239,0,0},{242,0,0},{244,0,0},{247,0,0},{250,0,0},{252,0,0},
    {255,0,0},{255,0,0},{255,3,0},{255,5,0},{255,8,0},{255,11,0},{255,13,0},
    {255,16,0},{255,19,0},{255,21,0},{255,24,0},{255,27,0},{255,30,0},
    {255,32,0},{255,35,0},{255,38,0},{255,40,0},{255,43,0},{255,46,0},
    {255,48,0},{255,51,0},{255,54,0},{255,56,0},{255,59,0},{255,62,0},
    {255,64,0},{255,67,0},{255,70,0},{255,72,0},{255,75,0},{255,78,0},
    {255,81,0},{255,83,0},{255,86,0},{255,89,0},{255,91,0},{255,94,0},
    {255,97,0},{255,99,0},{255,102,0},{255,105,0},{255,107,0},{255,110,0},
    {255,113,0},{255,115,0},{255,118,0},{255,121,0},{255,123,0},{255,126,0},
    {255,129,0},{255,132,0},{255,134,0},{255,137,0},{255,140,0},{255,142,0},
    {255,145,0},{255,148,0},{255,150,0},{255,153,0},{255,156,0},{255,158,0},
    {255,161,0},{255,164,0},{255,166,0},{255,169,0},{255,172,0},{255,174,0},
    {255,177,0},{255,180,0},{255,183,0},{255,185,0},{255,188,0},{255,191,0},
    {255,193,0},{255,196,0},{255,199,0},{255,201,0},{255,204,0},{255,207,0},
    {255,209,0},{255,212,0},{255,215,0},{255,217,0},{255,220,0},{255,223,0},
    {255,225,0},{255,228,0},{255,231,0},{255,234,0},{255,236,0},{255,239,0},
    {255,242,0},{255,244,0},{255,247,0},{255,250,0},{255,252,0},{255,255,0},
    {255,255,0},{255,255,4},{255,255,8},{255,255,12},{255,255,16},
    {255,255,20},{255,255,24},{255,255,28},{255,255,32},{255,255,36},
    {255,255,40},{255,255,45},{255,255,49},{255,255,53},{255,255,57},
    {255,255,61},{255,255,65},{255,255,69},{255,255,73},{255,255,77},
    {255,255,81},{255,255,85},{255,255,89},{255,255,93},{255,255,97},
    {255,255,101},{255,255,105},{255,255,109},{255,255,113},{255,255,117},
    {255,255,121},{255,255,125},{255,255,130},{255,255,134},{255,255,138},
    {255,255,142},{255,255,146},{255,255,150},{255,255,154},{255,255,158},
    {255,255,162},{255,255,166},{255,255,170},{255,255,174},{255,255,178},
    {255,255,182},{255,255,186},{255,255,190},{255,255,194},{255,255,198},
    {255,255,202},{255,255,206},{255,255,210},{255,255,215},{255,255,219},
    {255,255,223},{255,255,227},{255,255,231},{255,255,235},{255,255,239},
    {255,255,243},{255,255,247},{255,255,251},{255,255,255},
};

/* -----------------------------------------------------------------------
 * Write WAV file
 * -----------------------------------------------------------------------*/

static int write_wav(const char *path, const double *signal, unsigned n,
                     unsigned sample_rate)
{
    float *fdata = malloc(n * sizeof(float));
    if (!fdata) { fprintf(stderr, "allocation failed\n"); return -1; }
    for (unsigned i = 0; i < n; i++) fdata[i] = (float)signal[i];

    int ret = FIO_write_wav(path, fdata, n, sample_rate);
    free(fdata);
    if (ret != 0) fprintf(stderr, "Error writing %s\n", path);
    return ret;
}

/* -----------------------------------------------------------------------
 * Write interactive HTML spectrogram (Plotly heatmap)
 * -----------------------------------------------------------------------*/

static int write_html(const char *path, const double *spec_db,
                      unsigned num_frames, unsigned num_bins,
                      const double *times_s, const double *freqs_hz,
                      const char *text)
{
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); return -1; }

    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
             "Text: \"%s\" &nbsp;|&nbsp; 400–7300 Hz &nbsp;|&nbsp; "
             "2 s at 16 kHz", text);

    plot_html_begin(fp, "Spectrogram Text Art", subtitle, 0);

    fprintf(fp,
        "  <div id=\"spectrogram\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>How it works:</strong> Each character is rasterised with a\n"
        "    5&times;7 bitmap font. Each \"on\" pixel becomes a sine wave at the\n"
        "    corresponding frequency. The spectrogram reveals the hidden text.\n"
        "  </div>\n\n");

    fprintf(fp, "  <script>\n");
    plot_html_js_array(fp, "times", times_s, num_frames, "%.4f");
    plot_html_js_array(fp, "freqs", freqs_hz, num_bins, "%.2f");

    fprintf(fp, "    const z = [\n");
    for (unsigned k = 0; k < num_bins; k++) {
        fprintf(fp, "      [");
        for (unsigned f = 0; f < num_frames; f++) {
            fprintf(fp, "%.2f", spec_db[f * num_bins + k]);
            if (f + 1 < num_frames) fprintf(fp, ",");
        }
        fprintf(fp, "]");
        if (k + 1 < num_bins) fprintf(fp, ",");
        fprintf(fp, "\n");
    }
    fprintf(fp, "    ];\n\n");

    fprintf(fp,
        "    Plotly.newPlot('spectrogram', [{\n"
        "      type: 'heatmap',\n"
        "      x: times,\n"
        "      y: freqs,\n"
        "      z: z,\n"
        "      colorscale: 'Hot',\n"
        "      zmin: -100,\n"
        "      zmax: 0,\n"
        "      colorbar: { title: 'dB', thickness: 15 },\n"
        "      hovertemplate: 't: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>"
        "%%{z:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'Frequency (Hz)' },\n"
        "      margin: { t: 30, r: 80, b: 50, l: 70 },\n"
        "      height: 500\n"
        "    }, { responsive: true });\n"
        "  </script>\n");

    plot_html_end(fp);
    fclose(fp);
    return 0;
}

/* -----------------------------------------------------------------------
 * Write PNG spectrogram image
 * -----------------------------------------------------------------------*/

static int write_png(const char *path, const double *spec_db,
                     unsigned num_frames, unsigned num_bins,
                     int use_color)
{
    /* Image: width = num_frames, height = num_bins (flipped so high freq is top) */
    unsigned width  = num_frames;
    unsigned height = num_bins;
    unsigned char *pixels = malloc(width * height * 3);
    if (!pixels) { fprintf(stderr, "allocation failed\n"); return -1; }

    double db_min = -100.0;
    double db_max = 0.0;

    for (unsigned y = 0; y < height; y++) {
        /* Flip: row 0 in image = highest frequency bin */
        unsigned bin = height - 1 - y;
        for (unsigned x = 0; x < width; x++) {
            double db = spec_db[x * num_bins + bin];
            /* Clamp and normalise to [0, 1] */
            double t = (db - db_min) / (db_max - db_min);
            if (t < 0.0) t = 0.0;
            if (t > 1.0) t = 1.0;

            unsigned idx = (y * width + x) * 3;
            if (use_color) {
                unsigned ci = (unsigned)(t * 255.0);
                if (ci > 255) ci = 255;
                pixels[idx + 0] = hot_lut[ci][0];
                pixels[idx + 1] = hot_lut[ci][1];
                pixels[idx + 2] = hot_lut[ci][2];
            } else {
                unsigned char v = (unsigned char)(t * 255.0);
                pixels[idx + 0] = v;
                pixels[idx + 1] = v;
                pixels[idx + 2] = v;
            }
        }
    }

    int ret = stbi_write_png(path, (int)width, (int)height, 3, pixels,
                             (int)(width * 3));
    free(pixels);
    if (!ret) { fprintf(stderr, "Error writing %s\n", path); return -1; }
    return 0;
}

/* -----------------------------------------------------------------------
 * Main
 * -----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    /* Parse arguments */
    const char *text = "HELLO";
    int use_color = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--colormap") == 0 && i + 1 < argc) {
            if (strcmp(argv[i + 1], "grayscale") == 0)
                use_color = 0;
            i++;
        } else {
            text = argv[i];
        }
    }

    /* Synthesis parameters */
    const double sample_rate  = 16000.0;
    const double freq_lo      = 400.0;   /* 200 Hz margin above 0 Hz */
    const double freq_hi      = 7300.0;  /* 200 Hz margin below Nyquist region */
    const double duration_sec = 2.25;
    const double pad_sec      = 0.5;     /* silence before and after text */

    /* STFT parameters: large FFT, small hop for sharp spectrogram */
    const unsigned fft_n = 1024;
    const unsigned hop   = 16;

    /* Synthesise the audio into a temporary buffer, then pad with silence */
    unsigned max_text = (unsigned)(sample_rate * duration_sec) + 1024;
    double *text_buf = malloc(max_text * sizeof(double));
    if (!text_buf) { fprintf(stderr, "allocation failed\n"); return 1; }

    unsigned text_samples = MD_spectrogram_text(text_buf, max_text, text,
                                                freq_lo, freq_hi,
                                                duration_sec, sample_rate);

    unsigned pad_samples  = (unsigned)(pad_sec * sample_rate);
    unsigned num_samples  = pad_samples + text_samples + pad_samples;
    double *signal = calloc(num_samples, sizeof(double));
    if (!signal) { fprintf(stderr, "allocation failed\n"); return 1; }

    memcpy(signal + pad_samples, text_buf, text_samples * sizeof(double));
    free(text_buf);

    printf("Synthesised %u samples (%.3f s) for \"%s\" "
           "(%.0f ms silence padding)\n",
           num_samples, (double)num_samples / sample_rate, text,
           pad_sec * 1000.0);

    /* Write WAV */
    if (write_wav("spectrogram_text.wav", signal, num_samples,
                  (unsigned)sample_rate) != 0)
        return 1;
    printf("Wrote spectrogram_text.wav\n");

    /* Compute STFT */
    unsigned num_frames = MD_stft_num_frames(num_samples, fft_n, hop);
    unsigned num_bins   = fft_n / 2 + 1;

    double *mag     = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *spec_db = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *times_s = malloc(num_frames * sizeof(double));
    double *freqs   = malloc(num_bins * sizeof(double));

    if (!mag || !spec_db || !times_s || !freqs) {
        fprintf(stderr, "allocation failed\n"); return 1;
    }

    MD_stft(signal, num_samples, fft_n, hop, mag);

    /* Convert to dB */
    for (unsigned i = 0; i < num_frames * num_bins; i++)
        spec_db[i] = 20.0 * log10(fmax(mag[i] / (double)fft_n, 1e-6));

    /* Build axes */
    for (unsigned f = 0; f < num_frames; f++)
        times_s[f] = (double)(f * hop) / sample_rate;
    for (unsigned k = 0; k < num_bins; k++)
        freqs[k] = (double)k * sample_rate / (double)fft_n;

    /* Write HTML */
    if (write_html("spectrogram_text.html", spec_db, num_frames, num_bins,
                   times_s, freqs, text) == 0)
        printf("Wrote spectrogram_text.html\n");

    /* Write PNG */
    if (write_png("spectrogram_text.png", spec_db, num_frames, num_bins,
                  use_color) == 0)
        printf("Wrote spectrogram_text.png (%s)\n",
               use_color ? "hot" : "grayscale");

    printf("STFT: N=%u, hop=%u, %u frames, %u bins\n",
           fft_n, hop, num_frames, num_bins);

    /* Cleanup */
    free(freqs);
    free(times_s);
    free(spec_db);
    free(mag);
    free(signal);
    MD_shutdown();

    return 0;
}
