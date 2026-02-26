/**
 * @file spectrogram_text.c
 * @brief Example: generate audio whose spectrogram displays readable text.
 *
 * Usage:
 *   ./spectrogram_text [TEXT] [--colormap viridis|grayscale]
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
 * Viridis 256-entry colormap (R, G, B)
 * -----------------------------------------------------------------------*/

static const unsigned char viridis_lut[256][3] = {
    {68,1,84},{68,2,86},{69,4,87},{69,5,89},{70,7,90},{70,8,92},{70,10,93},
    {70,11,94},{71,13,96},{71,14,97},{71,16,99},{71,17,100},{71,19,101},
    {72,20,103},{72,22,104},{72,23,105},{72,24,106},{72,26,108},{72,27,109},
    {72,28,110},{72,29,111},{72,31,112},{72,32,113},{72,33,115},{72,35,116},
    {72,36,117},{72,37,118},{72,38,119},{72,40,120},{72,41,121},{71,42,122},
    {71,44,122},{71,45,123},{71,46,124},{71,47,125},{70,48,126},{70,50,126},
    {70,51,127},{69,52,128},{69,53,129},{69,55,129},{68,56,130},{68,57,131},
    {68,58,131},{67,60,132},{67,61,132},{66,62,133},{66,63,133},{66,64,134},
    {65,66,134},{65,67,135},{64,68,135},{64,69,136},{63,71,136},{63,72,137},
    {62,73,137},{62,74,137},{62,76,138},{61,77,138},{61,78,138},{60,79,139},
    {60,80,139},{59,82,139},{59,83,140},{58,84,140},{58,85,140},{57,86,140},
    {57,87,141},{56,88,141},{56,90,141},{55,91,141},{55,92,141},{54,93,142},
    {54,94,142},{53,95,142},{53,96,142},{52,97,142},{52,98,142},{51,99,142},
    {51,100,142},{50,101,142},{50,102,142},{49,103,142},{49,104,142},
    {49,105,142},{48,106,142},{48,107,142},{47,108,142},{47,109,142},
    {46,110,142},{46,111,142},{46,112,142},{45,113,142},{45,114,142},
    {44,115,142},{44,116,141},{44,117,141},{43,118,141},{43,119,141},
    {42,120,141},{42,121,141},{42,122,140},{41,123,140},{41,124,140},
    {41,125,140},{40,126,140},{40,127,139},{39,128,139},{39,129,139},
    {39,130,139},{38,131,138},{38,132,138},{37,133,138},{37,134,137},
    {37,135,137},{36,136,137},{36,137,136},{35,138,136},{35,139,136},
    {35,140,135},{34,141,135},{34,142,134},{34,143,134},{33,144,133},
    {33,145,133},{33,146,132},{32,147,132},{32,148,131},{32,149,131},
    {31,150,130},{31,151,130},{31,152,129},{31,153,128},{31,154,128},
    {30,155,127},{30,156,127},{30,157,126},{30,158,125},{30,159,124},
    {30,160,124},{30,161,123},{30,162,122},{31,163,122},{31,164,121},
    {31,165,120},{31,166,119},{32,167,119},{32,168,118},{33,169,117},
    {33,170,116},{34,171,115},{35,172,114},{35,173,114},{36,174,113},
    {37,175,112},{38,176,111},{39,177,110},{40,178,109},{41,179,108},
    {42,180,107},{44,181,106},{45,182,105},{46,183,104},{48,184,103},
    {49,185,102},{51,186,101},{52,187,100},{54,188,99},{56,189,97},
    {57,190,96},{59,191,95},{61,192,94},{63,193,92},{64,194,91},
    {66,195,90},{68,196,88},{70,197,87},{72,198,86},{74,199,84},
    {76,200,83},{78,201,81},{80,202,80},{82,203,78},{84,204,77},
    {86,205,75},{89,206,74},{91,207,72},{93,208,70},{95,209,69},
    {97,210,67},{100,211,65},{102,212,63},{105,213,62},{107,214,60},
    {109,215,58},{112,216,56},{114,217,54},{117,218,52},{119,219,50},
    {122,220,48},{124,221,46},{127,222,44},{129,222,42},{132,223,40},
    {134,224,38},{137,225,35},{139,226,33},{142,227,31},{145,228,29},
    {147,229,26},{150,230,24},{153,231,22},{155,232,19},{158,233,17},
    {161,234,14},{163,235,12},{166,235,10},{169,236,8},{171,237,6},
    {174,238,5},{177,239,4},{179,240,3},{182,241,3},{184,241,3},
    {187,242,3},{189,243,4},{192,244,5},{194,244,6},{197,245,8},
    {199,246,10},{202,247,13},{204,247,15},{207,248,18},{209,249,21},
    {212,250,24},{214,250,27},{217,251,30},{219,252,34},{222,253,37},
    {224,253,41},{226,254,44},{229,254,48},{231,255,52},{234,255,55},
    {236,255,59},{238,255,63},{241,255,67},{243,255,71},{245,255,75},
    {247,254,79},{249,254,83},{251,254,87},{253,254,91},
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
        "      colorscale: 'Viridis',\n"
        "      zmin: -80,\n"
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
                     int use_viridis)
{
    /* Image: width = num_frames, height = num_bins (flipped so high freq is top) */
    unsigned width  = num_frames;
    unsigned height = num_bins;
    unsigned char *pixels = malloc(width * height * 3);
    if (!pixels) { fprintf(stderr, "allocation failed\n"); return -1; }

    double db_min = -80.0;
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
            if (use_viridis) {
                unsigned ci = (unsigned)(t * 255.0);
                if (ci > 255) ci = 255;
                pixels[idx + 0] = viridis_lut[ci][0];
                pixels[idx + 1] = viridis_lut[ci][1];
                pixels[idx + 2] = viridis_lut[ci][2];
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
    int use_viridis = 1;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--colormap") == 0 && i + 1 < argc) {
            if (strcmp(argv[i + 1], "grayscale") == 0)
                use_viridis = 0;
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
                  use_viridis) == 0)
        printf("Wrote spectrogram_text.png (%s)\n",
               use_viridis ? "viridis" : "grayscale");

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
