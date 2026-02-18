/**
 * @file gen_signal_plots.c
 * @brief Docs-only utility: generate spectrum and spectrogram HTML plots
 *        for the signal generators guide.
 *
 * Generates 14 self-contained HTML files (7 magnitude spectra + 7 STFT
 * spectrograms) into guides/plots/ for embedding as iframes in the
 * Doxygen signal generators guide.
 *
 * It is NOT a user-facing example — it is invoked automatically by
 * `make docs`.
 *
 * Build (from repo root):
 *   make -C examples gen_signal_plots
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

/* Signal parameters: fs = N_SIGNAL = 8192 gives bin k = k Hz exactly,
 * zero spectral leakage, no windowing needed. */
#define SAMPLE_RATE  8192u
#define N_SIGNAL     8192u

/* STFT parameters: 256-point window (~31 ms at 8192 Hz), 75% overlap */
#define N_FFT   256u
#define HOP      64u

static void write_head(FILE *fp, const char *title)
{
    fprintf(fp,
        "<!DOCTYPE html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <title>%s</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n"
        "  <style>\n"
        "    * { box-sizing: border-box; margin: 0; padding: 0; }\n"
        "    body { font-family: system-ui, -apple-system, sans-serif;"
        " background: #fafafa; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <div id=\"plot\"></div>\n"
        "  <script>\n",
        title);
}

static void write_foot(FILE *fp)
{
    fprintf(fp,
        "  </script>\n"
        "</body>\n"
        "</html>\n");
}

/**
 * Write a magnitude spectrum HTML file.
 *
 * Applies the codebase-standard single-sided amplitude normalization:
 *   mag[k] /= N_SIGNAL
 *   if (k > 0 && k < N_SIGNAL/2) mag[k] *= 2.0
 * dB conversion is done client-side in JavaScript.
 */
static int write_spectrum_html(const char *path, const char *title,
                               const double *signal)
{
    const unsigned num_bins = N_SIGNAL / 2 + 1;
    double *mag = malloc(num_bins * sizeof(double));
    if (!mag) {
        fprintf(stderr, "allocation failed for %s\n", path);
        return -1;
    }

    MD_magnitude_spectrum(signal, N_SIGNAL, mag);

    /* Single-sided amplitude normalization */
    for (unsigned k = 0; k < num_bins; k++) {
        mag[k] /= (double)N_SIGNAL;
        if (k > 0 && k < N_SIGNAL / 2)
            mag[k] *= 2.0;
    }

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        free(mag);
        return -1;
    }

    write_head(fp, title);

    /* Frequency axis: bin k = k Hz exactly (fs = N_SIGNAL = 8192) */
    fprintf(fp, "    const freqs = Array.from({length: %u}, (_, k) => k);\n",
            num_bins);

    /* Embed magnitude values; dB conversion done in JS */
    plot_html_js_array(fp, "mags", mag, num_bins, "%.7g");

    fprintf(fp,
        "    const mags_db = mags.map(m => 20 * Math.log10(Math.max(m, 1e-6)));\n"
        "    Plotly.newPlot('plot', [{\n"
        "      x: freqs, y: mags_db,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#2563eb', width: 1.2 },\n"
        "      hovertemplate: '%%{x:.0f} Hz<br>%%{y:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Frequency (Hz)', range: [0, 4096] },\n"
        "      yaxis: { title: 'Magnitude (dB)', range: [-80, 5] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title);

    write_foot(fp);
    fclose(fp);
    free(mag);
    printf("  %s\n", path);
    return 0;
}

/**
 * Write a spectrogram HTML file using STFT (N_FFT=256, hop=64, 75% overlap).
 * dB conversion: 20*log10(mag / N_FFT), floored at -120 dB.
 */
static int write_spectrogram_html(const char *path, const char *title,
                                  const double *signal)
{
    const unsigned num_bins   = N_FFT / 2 + 1;
    const unsigned num_frames = MD_stft_num_frames(N_SIGNAL, N_FFT, HOP);

    double *mag_out = malloc((size_t)num_frames * num_bins * sizeof(double));
    if (!mag_out) {
        fprintf(stderr, "allocation failed for %s\n", path);
        return -1;
    }

    MD_stft(signal, N_SIGNAL, N_FFT, HOP, mag_out);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        free(mag_out);
        return -1;
    }

    write_head(fp, title);

    /* Time axis: frame f starts at f * HOP / SAMPLE_RATE seconds.
     * With HOP=64, SAMPLE_RATE=8192: step = 0.0078125 (exact) */
    const double time_step = (double)HOP / (double)SAMPLE_RATE;
    fprintf(fp,
            "    const times = Array.from({length: %u}, (_, f) => f * %.10g);\n",
            num_frames, time_step);

    /* Frequency axis: bin k = k * SAMPLE_RATE / N_FFT Hz.
     * With SAMPLE_RATE=8192, N_FFT=256: step = 32 Hz (exact integer) */
    const double freq_step = (double)SAMPLE_RATE / (double)N_FFT;
    fprintf(fp,
            "    const freqs = Array.from({length: %u}, (_, k) => k * %.10g);\n",
            num_bins, freq_step);

    /* Spectrogram matrix: z[k][f] — row index = freq bin (y-axis),
     * col index = time frame (x-axis), as required by Plotly heatmap. */
    fprintf(fp, "    const z = [\n");
    for (unsigned k = 0; k < num_bins; k++) {
        fprintf(fp, "      [");
        for (unsigned f = 0; f < num_frames; f++) {
            double db = 20.0 * log10(fmax(mag_out[f * num_bins + k] / (double)N_FFT, 1e-6));
            fprintf(fp, "%.1f", db);
            if (f + 1 < num_frames) fprintf(fp, ",");
        }
        fprintf(fp, "]");
        if (k + 1 < num_bins) fprintf(fp, ",");
        fprintf(fp, "\n");
    }
    fprintf(fp, "    ];\n");

    fprintf(fp,
        "    Plotly.newPlot('plot', [{\n"
        "      type: 'heatmap',\n"
        "      x: times, y: freqs, z: z,\n"
        "      colorscale: 'Viridis',\n"
        "      zmin: -80, zmax: 0,\n"
        "      colorbar: { title: 'dB', thickness: 12 },\n"
        "      hovertemplate: 't: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>%%{z:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'Frequency (Hz)' },\n"
        "      margin: { t: 35, r: 80, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title);

    write_foot(fp);
    fclose(fp);
    free(mag_out);
    printf("  %s\n", path);
    return 0;
}

int main(void)
{
    double *buf = malloc(N_SIGNAL * sizeof(double));
    if (!buf) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    printf("Generating signal plots (%u samples, %u Hz, 1.0s):\n",
           N_SIGNAL, SAMPLE_RATE);

    /* ----------------------------------------------------------------
     * Phase 1: magnitude spectra
     * All use the N=8192 FFT plan — created once on the first call,
     * reused for all 7 signals, released by MD_shutdown().
     * ----------------------------------------------------------------*/
    printf("Magnitude spectra:\n");

    MD_sine_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrum_html("guides/plots/sine_spectrum.html",
                        "Sine 440 Hz - Spectrum", buf);

    MD_square_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrum_html("guides/plots/square_spectrum.html",
                        "Square 440 Hz - Spectrum", buf);

    MD_sawtooth_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrum_html("guides/plots/sawtooth_spectrum.html",
                        "Sawtooth 440 Hz - Spectrum", buf);

    MD_chirp_linear(buf, N_SIGNAL, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_spectrum_html("guides/plots/chirp_linear_spectrum.html",
                        "Linear Chirp 20-4000 Hz - Spectrum", buf);

    MD_chirp_log(buf, N_SIGNAL, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_spectrum_html("guides/plots/chirp_log_spectrum.html",
                        "Log Chirp 20-4000 Hz - Spectrum", buf);

    MD_white_noise(buf, N_SIGNAL, 0.25, 42);
    write_spectrum_html("guides/plots/white_noise_spectrum.html",
                        "White Noise (sigma=0.25) - Spectrum", buf);

    memset(buf, 0, N_SIGNAL * sizeof(double));
    for (unsigned i = 0; i < 4; i++)
        buf[i * (N_SIGNAL / 4)] = 0.8;   /* positions 0, 2048, 4096, 6144 */
    write_spectrum_html("guides/plots/impulse_spectrum.html",
                        "Impulse Train (4 clicks) - Spectrum", buf);

    MD_shutdown();   /* release N=8192 plan before switching to N=256 */

    /* ----------------------------------------------------------------
     * Phase 2: spectrograms
     * All use the N=256 STFT plan — created once on the first call,
     * reused for all 7 signals, released by MD_shutdown().
     * ----------------------------------------------------------------*/
    printf("Spectrograms:\n");

    MD_sine_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/sine_spectrogram.html",
                           "Sine 440 Hz - Spectrogram", buf);

    MD_square_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/square_spectrogram.html",
                           "Square 440 Hz - Spectrogram", buf);

    MD_sawtooth_wave(buf, N_SIGNAL, 0.8, 440.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/sawtooth_spectrogram.html",
                           "Sawtooth 440 Hz - Spectrogram", buf);

    MD_chirp_linear(buf, N_SIGNAL, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/chirp_linear_spectrogram.html",
                           "Linear Chirp 20-4000 Hz - Spectrogram", buf);

    MD_chirp_log(buf, N_SIGNAL, 0.8, 20.0, 4000.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/chirp_log_spectrogram.html",
                           "Log Chirp 20-4000 Hz - Spectrogram", buf);

    MD_white_noise(buf, N_SIGNAL, 0.25, 42);
    write_spectrogram_html("guides/plots/white_noise_spectrogram.html",
                           "White Noise (sigma=0.25) - Spectrogram", buf);

    memset(buf, 0, N_SIGNAL * sizeof(double));
    for (unsigned i = 0; i < 4; i++)
        buf[i * (N_SIGNAL / 4)] = 0.8;
    write_spectrogram_html("guides/plots/impulse_spectrogram.html",
                           "Impulse Train (4 clicks) - Spectrogram", buf);

    MD_shutdown();   /* release N=256 STFT plan */

    free(buf);
    return 0;
}
