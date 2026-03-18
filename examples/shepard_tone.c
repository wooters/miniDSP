/**
 * @file shepard_tone.c
 * @brief Example: generate Shepard tone audio with spectrogram visualisation.
 *
 * Usage:
 *   ./shepard_tone [--rising | --falling | --static]
 *                  [--rate OCTAVES_PER_SEC]
 *                  [--octaves NUM]
 *                  [--base FREQ_HZ]
 *                  [--duration SEC]
 *
 * Outputs:
 *   shepard_tone.wav  — synthesised audio (WAV, IEEE float)
 *   shepard_tone.html — interactive Plotly spectrogram
 *
 * Build and run (from the repository root):
 *   make -C examples shepard_tone
 *   cd examples && ./shepard_tone
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "fileio.h"
#include "plot_html.h"

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
                      double rate, unsigned num_octaves, double base_freq)
{
    FILE *fp = fopen(path, "w");
    if (!fp) { fprintf(stderr, "cannot open %s\n", path); return -1; }

    const char *direction = (rate > 0) ? "Rising" :
                            (rate < 0) ? "Falling" : "Static";
    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
             "%s &nbsp;|&nbsp; base %.0f Hz &nbsp;|&nbsp; "
             "%u octaves &nbsp;|&nbsp; rate %.2f oct/s &nbsp;|&nbsp; 44.1 kHz",
             direction, base_freq, num_octaves, rate);

    plot_html_begin(fp, "Shepard Tone Spectrogram", subtitle, 0);

    fprintf(fp,
        "  <div id=\"spectrogram\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>Shepard tone:</strong> Multiple sine waves spaced one octave apart\n"
        "    glide continuously %s in pitch. A Gaussian spectral envelope\n"
        "    fades tones in at one end and out at the other, creating the auditory\n"
        "    illusion of endlessly %s pitch.\n"
        "  </div>\n\n",
        rate >= 0 ? "upward" : "downward",
        rate >= 0 ? "rising" : "falling");

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
        "      yaxis: { title: 'Frequency (Hz)', type: 'log' },\n"
        "      margin: { t: 30, r: 80, b: 50, l: 70 },\n"
        "      height: 500\n"
        "    }, { responsive: true });\n"
        "  </script>\n");

    plot_html_end(fp);
    fclose(fp);
    return 0;
}

/* -----------------------------------------------------------------------
 * Main
 * -----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
    MD_set_error_handler(NULL);  /* use default stderr handler */

    /* Defaults */
    double rate       = 0.5;    /* octaves per second (positive = rising) */
    unsigned num_oct  = 8;
    double base_freq  = 440.0;
    double duration   = 5.0;    /* seconds */
    double sample_rate = 44100.0;

    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--rising") == 0) {
            if (rate < 0) rate = -rate;
        } else if (strcmp(argv[i], "--falling") == 0) {
            if (rate > 0) rate = -rate;
        } else if (strcmp(argv[i], "--static") == 0) {
            rate = 0.0;
        } else if (strcmp(argv[i], "--rate") == 0 && i + 1 < argc) {
            rate = atof(argv[++i]);
        } else if (strcmp(argv[i], "--octaves") == 0 && i + 1 < argc) {
            num_oct = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--base") == 0 && i + 1 < argc) {
            base_freq = atof(argv[++i]);
        } else if (strcmp(argv[i], "--duration") == 0 && i + 1 < argc) {
            duration = atof(argv[++i]);
        }
    }

    unsigned N = (unsigned)(sample_rate * duration);

    printf("Shepard tone: %s, rate=%.2f oct/s, base=%.0f Hz, "
           "%u octaves, %.1f s, %.0f Hz\n",
           rate > 0 ? "rising" : rate < 0 ? "falling" : "static",
           rate, base_freq, num_oct, duration, sample_rate);

    /* Generate the Shepard tone */
    double *signal = malloc(N * sizeof(double));
    if (!signal) { fprintf(stderr, "allocation failed\n"); return 1; }

    /// [generate-signal]
    MD_shepard_tone(signal, N, 0.8, base_freq, sample_rate, rate, num_oct);
    /// [generate-signal]

    /* Write WAV */
    if (write_wav("shepard_tone.wav", signal, N, (unsigned)sample_rate) != 0)
        return 1;
    printf("Wrote shepard_tone.wav\n");

    /* Compute STFT for spectrogram */
    const unsigned fft_n = 2048;
    const unsigned hop   = 512;

    unsigned num_frames = MD_stft_num_frames(N, fft_n, hop);
    unsigned num_bins   = fft_n / 2 + 1;

    double *mag     = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *spec_db = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *times_s = malloc(num_frames * sizeof(double));
    double *freqs   = malloc(num_bins * sizeof(double));

    if (!mag || !spec_db || !times_s || !freqs) {
        fprintf(stderr, "allocation failed\n"); return 1;
    }

    MD_stft(signal, N, fft_n, hop, mag);

    /* Convert to dB */
    for (unsigned i = 0; i < num_frames * num_bins; i++)
        spec_db[i] = 20.0 * log10(fmax(mag[i] / (double)fft_n, 1e-6));

    /* Build axes */
    for (unsigned f = 0; f < num_frames; f++)
        times_s[f] = (double)(f * hop) / sample_rate;
    for (unsigned k = 0; k < num_bins; k++)
        freqs[k] = (double)k * sample_rate / (double)fft_n;

    /* Write HTML */
    if (write_html("shepard_tone.html", spec_db, num_frames, num_bins,
                   times_s, freqs, rate, num_oct, base_freq) == 0)
        printf("Wrote shepard_tone.html\n");

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
