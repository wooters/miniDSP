/**
 * @file gen_signal_plots.c
 * @brief Docs-only utility: generate HTML plots for docs guides.
 *
 * Generates:
 *   - 14 signal-generator plots (7 spectra + 7 spectrograms)
 *   - 6 simple-effects spectrograms (delay/tremolo/comb, before + after)
 *   - 8 window-function plots (4 time-domain + 4 spectra)
 *   - 8 FIR/convolution plots (4 responses + 4 spectra)
 *   - 3 pitch-detection plots (F0 tracks + ACF peak + FFT peak)
 *   - 5 mel/MFCC plots (input waveform + input spectrogram + filterbank
 *     shapes + mel energies + MFCC bars)
 *   - 1 DTMF spectrogram
 *   - 1 spectrogram-text spectrogram
 * into guides/plots/ for embedding as iframes in Doxygen guides.
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

/* Window visualisation parameters */
#define WINDOW_N        256u
#define WINDOW_FFT_VIS 4096u
#define FIR_TIME_SHOW    512u

/* Pitch-detection visualisation parameters */
#define PITCH_FRAME_N 1024u
#define PITCH_HOP      128u
#define PITCH_MIN_F0    80.0
#define PITCH_MAX_F0   400.0

/* Mel/MFCC visualisation parameters */
#define MEL_FRAME_N    1024u
#define MEL_NUM_MELS     26u
#define MEL_NUM_COEFFS   13u
#define MEL_MIN_FREQ    80.0
#define MEL_MAX_FREQ  3900.0

/* DTMF spectrogram parameters */
#define DTMF_SAMPLE_RATE 8000.0
#define DTMF_TONE_MS      70u
#define DTMF_PAUSE_MS     70u
#define DTMF_N_FFT       256u
#define DTMF_HOP           8u

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

/** Write a time-domain window plot: w[n] vs sample index n. */
static int write_window_time_html(const char *path, const char *title,
                                  const double *window, unsigned n)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);

    fprintf(fp, "    const idx = Array.from({length: %u}, (_, i) => i);\n", n);
    plot_html_js_array(fp, "w", window, n, "%.8g");

    fprintf(fp,
        "    Plotly.newPlot('plot', [{\n"
        "      x: idx, y: w,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#2563eb', width: 1.4 },\n"
        "      hovertemplate: 'n=%%{x}<br>w[n]=%%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Sample Index (n)', range: [0, %u] },\n"
        "      yaxis: { title: 'Amplitude', range: [-0.05, 1.05] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title, n - 1);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/**
 * Write a zero-padded window magnitude response:
 * one-sided amplitude in dB vs normalized frequency (cycles/sample).
 */
static int write_window_spectrum_html(const char *path, const char *title,
                                      const double *window, unsigned n_win,
                                      unsigned n_fft_vis)
{
    double *signal = calloc(n_fft_vis, sizeof(double));
    double *mag = malloc((n_fft_vis / 2 + 1) * sizeof(double));
    if (!signal || !mag) {
        fprintf(stderr, "allocation failed for %s\n", path);
        free(signal);
        free(mag);
        return -1;
    }

    memcpy(signal, window, n_win * sizeof(double));
    MD_magnitude_spectrum(signal, n_fft_vis, mag);

    /* Single-sided amplitude normalization */
    for (unsigned k = 0; k < n_fft_vis / 2 + 1; k++) {
        mag[k] /= (double)n_fft_vis;
        if (k > 0 && k < n_fft_vis / 2)
            mag[k] *= 2.0;
    }

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        free(signal);
        free(mag);
        return -1;
    }

    write_head(fp, title);

    fprintf(fp,
            "    const f = Array.from({length: %u}, (_, k) => k / %u.0);\n",
            n_fft_vis / 2 + 1, n_fft_vis);
    plot_html_js_array(fp, "mags", mag, n_fft_vis / 2 + 1, "%.8g");

    fprintf(fp,
        "    const mags_db = mags.map(m => 20 * Math.log10(Math.max(m, 1e-6)));\n"
        "    Plotly.newPlot('plot', [{\n"
        "      x: f, y: mags_db,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#2563eb', width: 1.2 },\n"
        "      hovertemplate: 'f: %%{x:.4f} cycles/sample<br>%%{y:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Normalized Frequency (cycles/sample)', range: [0, 0.5] },\n"
        "      yaxis: { title: 'Magnitude (dB)', range: [-120, 5] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title);

    write_foot(fp);
    fclose(fp);
    free(signal);
    free(mag);
    printf("  %s\n", path);
    return 0;
}

/** Write a time-domain signal plot (first show_n samples). */
static int write_signal_time_html(const char *path, const char *title,
                                  const double *signal, unsigned n,
                                  unsigned show_n)
{
    if (show_n > n) show_n = n;

    double ymin = signal[0];
    double ymax = signal[0];
    for (unsigned i = 1; i < show_n; i++) {
        if (signal[i] < ymin) ymin = signal[i];
        if (signal[i] > ymax) ymax = signal[i];
    }

    double yrange = ymax - ymin;
    double pad = (yrange > 1e-12) ? 0.10 * yrange : 0.1;
    double ylo = ymin - pad;
    double yhi = ymax + pad;

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);

    fprintf(fp, "    const idx = Array.from({length: %u}, (_, i) => i);\n", show_n);
    plot_html_js_array(fp, "sig", signal, show_n, "%.8g");

    fprintf(fp,
        "    Plotly.newPlot('plot', [{\n"
        "      x: idx, y: sig,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#059669', width: 1.3 },\n"
        "      hovertemplate: 'n=%%{x}<br>x[n]=%%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Sample Index (n)', range: [0, %u] },\n"
        "      yaxis: { title: 'Amplitude', range: [%.8g, %.8g] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title, show_n - 1, ylo, yhi);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/** Write an F0 track comparison plot (truth vs autocorrelation vs FFT). */
static int write_pitch_tracks_html(const char *path, const char *title,
                                   const double *times, const double *truth,
                                   const double *acf, const double *fft,
                                   unsigned num_frames)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);
    plot_html_js_array(fp, "times", times, num_frames, "%.8g");
    plot_html_js_array(fp, "f0Truth", truth, num_frames, "%.8g");
    plot_html_js_array(fp, "f0Acf", acf, num_frames, "%.8g");
    plot_html_js_array(fp, "f0Fft", fft, num_frames, "%.8g");

    fprintf(fp,
        "    Plotly.newPlot('plot', [\n"
        "      { x: times, y: f0Truth, name: 'Ground Truth', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#111827', width: 2.0 } },\n"
        "      { x: times, y: f0Acf, name: 'Autocorrelation F0', type: 'scatter', mode: 'lines+markers',\n"
        "        marker: { size: 4 }, line: { color: '#2563eb', width: 1.4 } },\n"
        "      { x: times, y: f0Fft, name: 'FFT F0', type: 'scatter', mode: 'lines+markers',\n"
        "        marker: { size: 4 }, line: { color: '#dc2626', width: 1.3 } }\n"
        "    ], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'F0 (Hz)' },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n",
        title);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/** Write an autocorrelation-vs-lag plot with selected F0 lag marker. */
static int write_pitch_acf_peak_html(const char *path, const char *title,
                                     const double *acf, unsigned lag_min,
                                     unsigned lag_max, double selected_lag)
{
    unsigned len = lag_max - lag_min + 1;
    double *lags = malloc(len * sizeof(double));
    double *vals = malloc(len * sizeof(double));
    if (!lags || !vals) {
        fprintf(stderr, "allocation failed for %s\n", path);
        free(vals);
        free(lags);
        return -1;
    }

    for (unsigned i = 0; i < len; i++) {
        lags[i] = (double)(lag_min + i);
        vals[i] = acf[lag_min + i];
    }

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        free(vals);
        free(lags);
        return -1;
    }

    write_head(fp, title);
    plot_html_js_array(fp, "lags", lags, len, "%.8g");
    plot_html_js_array(fp, "acfVals", vals, len, "%.8g");

    double selected_val = 0.0;
    if (selected_lag > 0.0) {
        unsigned idx = (unsigned)lround(selected_lag);
        if (idx >= lag_min && idx <= lag_max) {
            selected_val = acf[idx];
        }
    }

    fprintf(fp,
        "    Plotly.newPlot('plot', [\n"
        "      { x: lags, y: acfVals, name: 'Autocorrelation', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#2563eb', width: 1.4 } },\n"
        "      { x: [%.8g], y: [%.8g], name: 'Selected Lag', mode: 'markers',\n"
        "        marker: { color: '#dc2626', size: 9, symbol: 'diamond' } }\n"
        "    ], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Lag (samples)' },\n"
        "      yaxis: { title: 'Normalised R[lag]', range: [-1.1, 1.1] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n",
        selected_lag, selected_val, title);

    write_foot(fp);
    fclose(fp);
    free(vals);
    free(lags);
    printf("  %s\n", path);
    return 0;
}

/** Write an FFT magnitude plot with selected peak marker. */
static int write_pitch_fft_peak_html(const char *path, const char *title,
                                     const double *mag, unsigned num_bins,
                                     double selected_freq_hz)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);
    plot_html_js_array(fp, "magVals", mag, num_bins, "%.8g");

    const double freq_step = (double)SAMPLE_RATE / (double)PITCH_FRAME_N;

    double selected_mag_db = -120.0;
    if (selected_freq_hz > 0.0) {
        unsigned k = (unsigned)lround(selected_freq_hz / freq_step);
        if (k < num_bins) {
            selected_mag_db = 20.0 * log10(fmax(mag[k], 1e-6));
        }
    }

    fprintf(fp,
        "    const freqs = Array.from({length: %u}, (_, k) => k * %.10g);\n"
        "    const magsDb = magVals.map(m => 20 * Math.log10(Math.max(m, 1e-6)));\n"
        "    Plotly.newPlot('plot', [\n"
        "      { x: freqs, y: magsDb, name: 'Magnitude (dB)', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#dc2626', width: 1.3 } },\n"
        "      { x: [%.8g], y: [%.8g], name: 'Selected Peak', mode: 'markers',\n"
        "        marker: { color: '#2563eb', size: 9, symbol: 'diamond' } }\n"
        "    ], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Frequency (Hz)', range: [0, %u] },\n"
        "      yaxis: { title: 'Magnitude (dB)', range: [-100, 10] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n",
        num_bins, freq_step, selected_freq_hz, selected_mag_db, title, SAMPLE_RATE / 2);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/** Write mel filterbank triangular shapes over frequency. */
static int write_mel_filterbank_html(const char *path, const char *title,
                                     const double *fb, unsigned num_mels,
                                     unsigned N, double sample_rate)
{
    unsigned num_bins = N / 2 + 1;
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);
    const double freq_step = sample_rate / (double)N;

    fprintf(fp,
            "    const freqs = Array.from({length: %u}, (_, k) => k * %.10g);\n"
            "    const traces = [];\n",
            num_bins, freq_step);

    for (unsigned m = 0; m < num_mels; m++) {
        const double *row = fb + (size_t)m * num_bins;
        fprintf(fp, "    traces.push({ x: freqs, y: [");
        for (unsigned k = 0; k < num_bins; k++) {
            fprintf(fp, "%.7g", row[k]);
            if (k + 1 < num_bins) fprintf(fp, ",");
        }
        fprintf(fp,
            "], type: 'scatter', mode: 'lines',\n"
            "      line: { width: 1.2 }, opacity: 0.85,\n"
            "      showlegend: false,\n"
            "      hovertemplate: 'Band %u<br>f: %%{x:.1f} Hz<br>w: %%{y:.3f}<extra></extra>'\n"
            "    });\n", m);
    }

    fprintf(fp,
        "    Plotly.newPlot('plot', traces, {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Frequency (Hz)', range: [0, %.8g] },\n"
        "      yaxis: { title: 'Weight', range: [0, 1.05] },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title, sample_rate * 0.5);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/** Write a mel-energy bar chart for one analysis frame. */
static int write_mel_energies_html(const char *path, const char *title,
                                   const double *centers_hz,
                                   const double *mel, unsigned num_mels)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);
    plot_html_js_array(fp, "centers", centers_hz, num_mels, "%.8g");
    plot_html_js_array(fp, "melVals", mel, num_mels, "%.10g");

    fprintf(fp,
        "    Plotly.newPlot('plot', [{\n"
        "      x: centers, y: melVals,\n"
        "      type: 'bar', marker: { color: '#2563eb' },\n"
        "      hovertemplate: 'Center: %%{x:.1f} Hz<br>Energy: %%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Approximate Mel Band Center (Hz)' },\n"
        "      yaxis: { title: 'Mel Energy', rangemode: 'tozero' },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 70 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        title);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

/** Write a bar chart of MFCC coefficients (C0 included). */
static int write_mfcc_html(const char *path, const char *title,
                           const double *mfcc, unsigned num_coeffs)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    write_head(fp, title);
    plot_html_js_array(fp, "coeffs", mfcc, num_coeffs, "%.10g");

    fprintf(fp,
        "    const idx = Array.from({length: %u}, (_, i) => i);\n"
        "    Plotly.newPlot('plot', [{\n"
        "      x: idx, y: coeffs,\n"
        "      type: 'bar', marker: { color: '#dc2626' },\n"
        "      hovertemplate: 'C%%{x}: %%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: '%s', font: { size: 13 } },\n"
        "      xaxis: { title: 'Coefficient Index' },\n"
        "      yaxis: { title: 'Value' },\n"
        "      margin: { t: 35, r: 20, b: 50, l: 60 },\n"
        "      height: 360\n"
        "    }, { responsive: true });\n",
        num_coeffs, title);

    write_foot(fp);
    fclose(fp);
    printf("  %s\n", path);
    return 0;
}

int main(void)
{
    double *buf = malloc(N_SIGNAL * sizeof(double));
    double *fx = malloc(N_SIGNAL * sizeof(double));
    if (!buf || !fx) {
        fprintf(stderr, "allocation failed\n");
        free(fx);
        free(buf);
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

    /* ----------------------------------------------------------------
     * Phase 2b: simple effects spectrograms (before/after)
     * ----------------------------------------------------------------*/
    printf("Simple effect spectrograms:\n");

    /* Delay/echo: click train source */
    memset(buf, 0, N_SIGNAL * sizeof(double));
    const unsigned click_step = SAMPLE_RATE * 35u / 100u;  /* 0.35 s */
    for (unsigned i = 0; i < N_SIGNAL; i += click_step) {
        buf[i] = 0.9;
    }
    write_spectrogram_html("guides/plots/effect_delay_before_spectrogram.html",
                           "Delay/Echo - Before (dry clicks) Spectrogram", buf);

    MD_delay_echo(buf, fx, N_SIGNAL, SAMPLE_RATE / 4u, 0.45, 1.0, 0.6); /* 250 ms */
    write_spectrogram_html("guides/plots/effect_delay_after_spectrogram.html",
                           "Delay/Echo - After Spectrogram", fx);

    /* Tremolo: steady sine tone source */
    MD_sine_wave(buf, N_SIGNAL, 0.8, 220.0, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/effect_tremolo_before_spectrogram.html",
                           "Tremolo - Before (220 Hz sine) Spectrogram", buf);

    MD_tremolo(buf, fx, N_SIGNAL, 5.0, 0.8, SAMPLE_RATE);
    write_spectrogram_html("guides/plots/effect_tremolo_after_spectrogram.html",
                           "Tremolo - After Spectrogram", fx);

    /* Comb reverb: short decaying tone burst */
    memset(buf, 0, N_SIGNAL * sizeof(double));
    const unsigned burst_len = SAMPLE_RATE / 5u;  /* 200 ms */
    for (unsigned i = 0; i < burst_len; i++) {
        double env = exp(-6.0 * (double)i / (double)burst_len);
        buf[i] = 0.85 * env * sin(2.0 * M_PI * 330.0 * (double)i / (double)SAMPLE_RATE);
    }
    write_spectrogram_html("guides/plots/effect_comb_before_spectrogram.html",
                           "Comb Reverb - Before (decaying burst) Spectrogram", buf);

    MD_comb_reverb(buf, fx, N_SIGNAL, SAMPLE_RATE * 3u / 100u, 0.75, 0.7, 0.6); /* 30 ms */
    write_spectrogram_html("guides/plots/effect_comb_after_spectrogram.html",
                           "Comb Reverb - After Spectrogram", fx);

    MD_shutdown();   /* release N=256 STFT plan */

    /* ----------------------------------------------------------------
     * Phase 3: window function visuals
     * All spectra use N=4096 zero-padded FFT for consistent comparison.
     * ----------------------------------------------------------------*/
    printf("Window function visuals:\n");

    double hann[WINDOW_N], hamming[WINDOW_N], blackman[WINDOW_N], rect[WINDOW_N];
    MD_Gen_Hann_Win(hann, WINDOW_N);
    MD_Gen_Hamming_Win(hamming, WINDOW_N);
    MD_Gen_Blackman_Win(blackman, WINDOW_N);
    MD_Gen_Rect_Win(rect, WINDOW_N);

    write_window_time_html("guides/plots/hann_window_time.html",
                           "Hanning Window - Time Domain", hann, WINDOW_N);
    write_window_spectrum_html("guides/plots/hann_window_spectrum.html",
                               "Hanning Window - Magnitude Response", hann,
                               WINDOW_N, WINDOW_FFT_VIS);

    write_window_time_html("guides/plots/hamming_window_time.html",
                           "Hamming Window - Time Domain", hamming, WINDOW_N);
    write_window_spectrum_html("guides/plots/hamming_window_spectrum.html",
                               "Hamming Window - Magnitude Response", hamming,
                               WINDOW_N, WINDOW_FFT_VIS);

    write_window_time_html("guides/plots/blackman_window_time.html",
                           "Blackman Window - Time Domain", blackman, WINDOW_N);
    write_window_spectrum_html("guides/plots/blackman_window_spectrum.html",
                               "Blackman Window - Magnitude Response", blackman,
                               WINDOW_N, WINDOW_FFT_VIS);

    write_window_time_html("guides/plots/rect_window_time.html",
                           "Rectangular Window - Time Domain", rect, WINDOW_N);
    write_window_spectrum_html("guides/plots/rect_window_spectrum.html",
                               "Rectangular Window - Magnitude Response", rect,
                               WINDOW_N, WINDOW_FFT_VIS);

    MD_shutdown();   /* release N=4096 plan */

    /* ----------------------------------------------------------------
     * Phase 4: FIR / convolution visuals
     * ----------------------------------------------------------------*/
    printf("FIR/convolution visuals:\n");

    double *work_a = malloc(N_SIGNAL * sizeof(double));
    double *work_b = malloc(N_SIGNAL * sizeof(double));
    if (!work_a || !work_b) {
        fprintf(stderr, "allocation failed\n");
        free(work_b);
        free(work_a);
        free(buf);
        return 1;
    }

    /* 4a) Time-domain direct convolution (full output length = N_SIGNAL) */
    const unsigned conv_kernel_len = 64;
    const unsigned conv_signal_len = N_SIGNAL - conv_kernel_len + 1;
    double *conv_in = malloc(conv_signal_len * sizeof(double));
    double conv_kernel[conv_kernel_len];
    if (!conv_in) {
        fprintf(stderr, "allocation failed\n");
        free(work_b);
        free(work_a);
        free(buf);
        return 1;
    }

    memset(conv_in, 0, conv_signal_len * sizeof(double));
    for (unsigned i = 128; i < conv_signal_len; i += 512) {
        conv_in[i] = 1.0;
    }
    for (unsigned k = 0; k < conv_kernel_len; k++) {
        conv_kernel[k] = exp(-0.08 * (double)k);
    }

    MD_convolution_time(conv_in, conv_signal_len, conv_kernel, conv_kernel_len, work_a);
    write_signal_time_html("guides/plots/conv_time_response.html",
                           "Time-Domain Convolution - Response", work_a,
                           N_SIGNAL, FIR_TIME_SHOW);
    write_spectrum_html("guides/plots/conv_time_spectrum.html",
                        "Time-Domain Convolution - Spectrum", work_a);

    /* 4b) Moving-average filter */
    MD_square_wave(buf, N_SIGNAL, 0.8, 220.0, SAMPLE_RATE);
    MD_white_noise(work_b, N_SIGNAL, 0.08, 123);
    for (unsigned i = 0; i < N_SIGNAL; i++) {
        buf[i] += work_b[i];
    }
    MD_moving_average(buf, N_SIGNAL, 16, work_a);
    write_signal_time_html("guides/plots/moving_average_response.html",
                           "Moving-Average Filter - Response", work_a,
                           N_SIGNAL, FIR_TIME_SHOW);
    write_spectrum_html("guides/plots/moving_average_spectrum.html",
                        "Moving-Average Filter - Spectrum", work_a);

    /* 4c) General FIR filter */
    const double fir_coeffs[] = {0.05, 0.09, 0.12, 0.15, 0.18, 0.15, 0.12, 0.09, 0.05};
    MD_fir_filter(buf, N_SIGNAL, fir_coeffs, 9, work_a);
    write_signal_time_html("guides/plots/fir_general_response.html",
                           "General FIR Filter - Response", work_a,
                           N_SIGNAL, FIR_TIME_SHOW);
    write_spectrum_html("guides/plots/fir_general_spectrum.html",
                        "General FIR Filter - Spectrum", work_a);

    /* 4d) FFT overlap-add fast convolution (same output as 4a kernel/input) */
    MD_convolution_fft_ola(conv_in, conv_signal_len, conv_kernel, conv_kernel_len, work_a);
    write_signal_time_html("guides/plots/conv_fft_ola_response.html",
                           "FFT Overlap-Add Convolution - Response", work_a,
                           N_SIGNAL, FIR_TIME_SHOW);
    write_spectrum_html("guides/plots/conv_fft_ola_spectrum.html",
                        "FFT Overlap-Add Convolution - Spectrum", work_a);

    /* ----------------------------------------------------------------
     * Phase 5: pitch-detection visuals
     * ----------------------------------------------------------------*/
    printf("Pitch detection visuals:\n");

    /* Build a voiced signal with piecewise F0 segments. */
    MD_white_noise(work_b, N_SIGNAL, 0.08, 777);
    const unsigned seg1 = N_SIGNAL / 3;
    const unsigned seg2 = 2 * N_SIGNAL / 3;
    double phase = 0.0;
    for (unsigned n = 0; n < N_SIGNAL; n++) {
        double f0 = (n < seg1) ? 140.0 : (n < seg2 ? 220.0 : 320.0);
        phase += 2.0 * M_PI * f0 / (double)SAMPLE_RATE;
        buf[n] = 0.75 * sin(phase)
               + 0.22 * sin(2.0 * phase)
               + 0.12 * sin(3.0 * phase)
               + work_b[n];
    }

    unsigned pitch_frames = MD_stft_num_frames(N_SIGNAL, PITCH_FRAME_N, PITCH_HOP);
    double *pitch_t = malloc(pitch_frames * sizeof(double));
    double *pitch_truth = malloc(pitch_frames * sizeof(double));
    double *pitch_acf = malloc(pitch_frames * sizeof(double));
    double *pitch_fft = malloc(pitch_frames * sizeof(double));
    if (!pitch_t || !pitch_truth || !pitch_acf || !pitch_fft) {
        fprintf(stderr, "allocation failed\n");
        free(pitch_fft);
        free(pitch_acf);
        free(pitch_truth);
        free(pitch_t);
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }

    for (unsigned f = 0; f < pitch_frames; f++) {
        unsigned start = f * PITCH_HOP;
        const double *frame = buf + start;
        unsigned center = start + PITCH_FRAME_N / 2;
        if (center >= N_SIGNAL) center = N_SIGNAL - 1;

        pitch_t[f] = (double)center / (double)SAMPLE_RATE;
        pitch_truth[f] = (center < seg1) ? 140.0 : (center < seg2 ? 220.0 : 320.0);
        pitch_acf[f] = MD_f0_autocorrelation(frame, PITCH_FRAME_N, SAMPLE_RATE,
                                             PITCH_MIN_F0, PITCH_MAX_F0);
        pitch_fft[f] = MD_f0_fft(frame, PITCH_FRAME_N, SAMPLE_RATE,
                                 PITCH_MIN_F0, PITCH_MAX_F0);
    }

    write_pitch_tracks_html("guides/plots/pitch_f0_tracks.html",
                            "Pitch Detection - Ground Truth vs Estimated F0",
                            pitch_t, pitch_truth, pitch_acf, pitch_fft, pitch_frames);

    /* Visualise one representative frame (middle segment around 220 Hz). */
    unsigned frame_idx = pitch_frames / 2;
    const double *frame = buf + frame_idx * PITCH_HOP;

    unsigned lag_min = (unsigned)floor((double)SAMPLE_RATE / PITCH_MAX_F0);
    unsigned lag_max = (unsigned)ceil((double)SAMPLE_RATE / PITCH_MIN_F0);
    if (lag_min < 1) lag_min = 1;
    if (lag_max > PITCH_FRAME_N - 1) lag_max = PITCH_FRAME_N - 1;

    double *acf = malloc((lag_max + 1) * sizeof(double));
    if (!acf) {
        fprintf(stderr, "allocation failed\n");
        free(pitch_fft);
        free(pitch_acf);
        free(pitch_truth);
        free(pitch_t);
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }
    MD_autocorrelation(frame, PITCH_FRAME_N, acf, lag_max + 1);
    double f0_acf_one = MD_f0_autocorrelation(frame, PITCH_FRAME_N, SAMPLE_RATE,
                                              PITCH_MIN_F0, PITCH_MAX_F0);
    double selected_lag = (f0_acf_one > 0.0) ? ((double)SAMPLE_RATE / f0_acf_one) : 0.0;
    write_pitch_acf_peak_html("guides/plots/pitch_acf_peak_frame.html",
                              "Pitch Detection - Autocorrelation Peak (One Frame)",
                              acf, lag_min, lag_max, selected_lag);
    free(acf);

    unsigned num_bins = PITCH_FRAME_N / 2 + 1;
    double *window = malloc(PITCH_FRAME_N * sizeof(double));
    double *frame_win = malloc(PITCH_FRAME_N * sizeof(double));
    double *mag = malloc(num_bins * sizeof(double));
    if (!window || !frame_win || !mag) {
        fprintf(stderr, "allocation failed\n");
        free(mag);
        free(frame_win);
        free(window);
        free(pitch_fft);
        free(pitch_acf);
        free(pitch_truth);
        free(pitch_t);
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }

    MD_Gen_Hann_Win(window, PITCH_FRAME_N);
    for (unsigned n = 0; n < PITCH_FRAME_N; n++) {
        frame_win[n] = frame[n] * window[n];
    }
    MD_magnitude_spectrum(frame_win, PITCH_FRAME_N, mag);
    for (unsigned k = 0; k < num_bins; k++) {
        mag[k] /= (double)PITCH_FRAME_N;
        if (k > 0 && k < PITCH_FRAME_N / 2) mag[k] *= 2.0;
    }
    double f0_fft_one = MD_f0_fft(frame, PITCH_FRAME_N, SAMPLE_RATE,
                                  PITCH_MIN_F0, PITCH_MAX_F0);
    write_pitch_fft_peak_html("guides/plots/pitch_fft_peak_frame.html",
                              "Pitch Detection - FFT Peak Pick (One Frame)",
                              mag, num_bins, f0_fft_one);
    free(mag);
    free(frame_win);
    free(window);

    free(pitch_fft);
    free(pitch_acf);
    free(pitch_truth);
    free(pitch_t);

    /* ----------------------------------------------------------------
     * Phase 6: mel/MFCC visuals
     * ----------------------------------------------------------------*/
    printf("Mel/MFCC visuals:\n");

    double mel_frame[MEL_FRAME_N];
    double mel_vals[MEL_NUM_MELS];
    double mfcc_vals[MEL_NUM_COEFFS];
    unsigned mel_bins = MEL_FRAME_N / 2 + 1;
    double *mel_fb = malloc((size_t)MEL_NUM_MELS * mel_bins * sizeof(double));
    double *mel_centers = malloc(MEL_NUM_MELS * sizeof(double));
    if (!mel_fb || !mel_centers) {
        fprintf(stderr, "allocation failed\n");
        free(mel_centers);
        free(mel_fb);
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }

    /* Deterministic signal used by all mel/MFCC visualisations. */
    for (unsigned n = 0; n < N_SIGNAL; n++) {
        double t = (double)n / (double)SAMPLE_RATE;
        buf[n] = 0.7 * sin(2.0 * M_PI * 440.0 * t)
               + 0.2 * cos(2.0 * M_PI * 1000.0 * t)
               + 0.1 * sin(2.0 * M_PI * 3000.0 * t);
    }
    write_signal_time_html("guides/plots/mel_input_waveform.html",
                           "Mel/MFCC Input Signal - Waveform",
                           buf, N_SIGNAL, 512);
    write_spectrogram_html("guides/plots/mel_input_spectrogram.html",
                           "Mel/MFCC Input Signal - Spectrogram",
                           buf);

    /* Use the first analysis frame from the deterministic signal. */
    memcpy(mel_frame, buf, MEL_FRAME_N * sizeof(double));

    MD_mel_filterbank(MEL_FRAME_N, (double)SAMPLE_RATE, MEL_NUM_MELS,
                      MEL_MIN_FREQ, MEL_MAX_FREQ, mel_fb);
    MD_mel_energies(mel_frame, MEL_FRAME_N, (double)SAMPLE_RATE, MEL_NUM_MELS,
                    MEL_MIN_FREQ, MEL_MAX_FREQ, mel_vals);
    MD_mfcc(mel_frame, MEL_FRAME_N, (double)SAMPLE_RATE,
            MEL_NUM_MELS, MEL_NUM_COEFFS, MEL_MIN_FREQ, MEL_MAX_FREQ, mfcc_vals);

    for (unsigned m = 0; m < MEL_NUM_MELS; m++) {
        const double *row = mel_fb + (size_t)m * mel_bins;
        double sum_w = 0.0;
        double sum_fw = 0.0;
        for (unsigned k = 0; k < mel_bins; k++) {
            double f_hz = (double)k * (double)SAMPLE_RATE / (double)MEL_FRAME_N;
            sum_w += row[k];
            sum_fw += f_hz * row[k];
        }
        mel_centers[m] = (sum_w > 0.0) ? (sum_fw / sum_w) : 0.0;
    }

    write_mel_filterbank_html("guides/plots/mel_filterbank_shapes.html",
                              "Mel Filterbank Shapes (HTK Mapping)",
                              mel_fb, MEL_NUM_MELS, MEL_FRAME_N, (double)SAMPLE_RATE);
    write_mel_energies_html("guides/plots/mel_energies_frame.html",
                            "Mel Energies (One Frame Example)",
                            mel_centers, mel_vals, MEL_NUM_MELS);
    write_mfcc_html("guides/plots/mfcc_frame.html",
                    "MFCCs (One Frame, C0 Included)",
                    mfcc_vals, MEL_NUM_COEFFS);

    free(mel_centers);
    free(mel_fb);

    MD_shutdown();   /* release mel/MFCC plans before DTMF phase */

    /* ----------------------------------------------------------------
     * Phase 7: DTMF spectrogram
     * Uses 8000 Hz sample rate (standard telephony) with N=256, hop=8
     * for high time-resolution visualization.
     * ----------------------------------------------------------------*/
    printf("DTMF spectrogram:\n");

    const char *dtmf_digits = "159#";
    unsigned num_dtmf = (unsigned)strlen(dtmf_digits);
    unsigned dtmf_len = MD_dtmf_signal_length(num_dtmf, DTMF_SAMPLE_RATE,
                                              DTMF_TONE_MS, DTMF_PAUSE_MS);
    double *dtmf_sig = malloc(dtmf_len * sizeof(double));
    if (!dtmf_sig) {
        fprintf(stderr, "allocation failed for DTMF signal\n");
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }
    MD_dtmf_generate(dtmf_sig, dtmf_digits, DTMF_SAMPLE_RATE,
                     DTMF_TONE_MS, DTMF_PAUSE_MS);

    /* Apply a 10 ms raised-cosine fade to each tone edge for a cleaner
     * spectrogram.  This is purely cosmetic — the library generates
     * standard rectangular-envelope DTMF tones. */
    {
        unsigned tone_samp  = (unsigned)(DTMF_TONE_MS  * DTMF_SAMPLE_RATE / 1000.0);
        unsigned pause_samp = (unsigned)(DTMF_PAUSE_MS * DTMF_SAMPLE_RATE / 1000.0);
        unsigned ramp = (unsigned)(0.010 * DTMF_SAMPLE_RATE);
        if (ramp > tone_samp / 2) ramp = tone_samp / 2;
        unsigned off = 0;
        for (unsigned d = 0; d < num_dtmf; d++) {
            for (unsigned i = 0; i < ramp; i++) {
                double g = 0.5 * (1.0 - cos(M_PI * i / ramp));
                dtmf_sig[off + i] *= g;
                dtmf_sig[off + tone_samp - 1 - i] *= g;
            }
            off += tone_samp + pause_samp;
        }
    }

    const unsigned dtmf_bins   = DTMF_N_FFT / 2 + 1;
    const unsigned dtmf_frames = MD_stft_num_frames(dtmf_len, DTMF_N_FFT, DTMF_HOP);

    double *dtmf_mag = malloc((size_t)dtmf_frames * dtmf_bins * sizeof(double));
    if (!dtmf_mag) {
        fprintf(stderr, "allocation failed for DTMF STFT\n");
        free(dtmf_sig);
        free(conv_in);
        free(work_b);
        free(work_a);
        free(fx);
        free(buf);
        return 1;
    }
    MD_stft(dtmf_sig, dtmf_len, DTMF_N_FFT, DTMF_HOP, dtmf_mag);

    {
        const char *path  = "guides/plots/dtmf_spectrogram.html";
        const char *title = "DTMF Sequence \"159#\" - Spectrogram";
        FILE *fp = fopen(path, "w");
        if (!fp) {
            fprintf(stderr, "cannot open %s for writing\n", path);
            free(dtmf_mag);
            free(dtmf_sig);
            free(conv_in);
            free(work_b);
            free(work_a);
            free(fx);
            free(buf);
            return 1;
        }

        write_head(fp, title);

        /* Time axis */
        const double dtmf_time_step = (double)DTMF_HOP / DTMF_SAMPLE_RATE;
        fprintf(fp,
                "    const times = Array.from({length: %u}, (_, f) => f * %.10g);\n",
                dtmf_frames, dtmf_time_step);

        /* Frequency axis */
        const double dtmf_freq_step = DTMF_SAMPLE_RATE / (double)DTMF_N_FFT;
        fprintf(fp,
                "    const freqs = Array.from({length: %u}, (_, k) => k * %.10g);\n",
                dtmf_bins, dtmf_freq_step);

        /* Spectrogram matrix: z[k][f] */
        fprintf(fp, "    const z = [\n");
        for (unsigned k = 0; k < dtmf_bins; k++) {
            fprintf(fp, "      [");
            for (unsigned f = 0; f < dtmf_frames; f++) {
                double db = 20.0 * log10(fmax(dtmf_mag[f * dtmf_bins + k]
                                              / (double)DTMF_N_FFT, 1e-6));
                fprintf(fp, "%.1f", db);
                if (f + 1 < dtmf_frames) fprintf(fp, ",");
            }
            fprintf(fp, "]");
            if (k + 1 < dtmf_bins) fprintf(fp, ",");
            fprintf(fp, "\n");
        }
        fprintf(fp, "    ];\n");

        /* DTMF reference frequencies for horizontal lines */
        fprintf(fp,
            "    const dtmfFreqs = [\n"
            "      {f: 697,  label: '697 Hz'},  {f: 770,  label: '770 Hz'},\n"
            "      {f: 852,  label: '852 Hz'},  {f: 941,  label: '941 Hz'},\n"
            "      {f: 1209, label: '1209 Hz'}, {f: 1336, label: '1336 Hz'},\n"
            "      {f: 1477, label: '1477 Hz'}, {f: 1633, label: '1633 Hz'}\n"
            "    ];\n");

        /* Build Plotly shapes array for reference lines */
        fprintf(fp,
            "    const shapes = dtmfFreqs.map(d => ({\n"
            "      type: 'line', xref: 'paper', x0: 0, x1: 1,\n"
            "      yref: 'y', y0: d.f, y1: d.f,\n"
            "      line: { color: 'rgba(255,255,255,0.5)', width: 1, dash: 'dot' }\n"
            "    }));\n");

        /* Build annotations for frequency labels */
        fprintf(fp,
            "    const annotations = dtmfFreqs.map(d => ({\n"
            "      xref: 'paper', x: 1.01, yref: 'y', y: d.f,\n"
            "      text: d.label, showarrow: false,\n"
            "      font: { size: 9, color: '#666' }, xanchor: 'left'\n"
            "    }));\n");

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
            "      yaxis: { title: 'Frequency (Hz)', range: [0, 2000] },\n"
            "      margin: { t: 35, r: 80, b: 50, l: 60 },\n"
            "      height: 360,\n"
            "      shapes: shapes,\n"
            "      annotations: annotations\n"
            "    }, { responsive: true });\n",
            title);

        write_foot(fp);
        fclose(fp);
        printf("  %s\n", path);
    }

    free(dtmf_mag);
    free(dtmf_sig);

    MD_shutdown();   /* release DTMF STFT plan */

    /* ----------------------------------------------------------------
     * Phase 8: Spectrogram text spectrogram
     * Uses 16000 Hz sample rate, N=1024, hop=16 to match the example
     * program's documented recipe.
     * ----------------------------------------------------------------*/
#define SPECTEXT_SAMPLE_RATE 16000.0
#define SPECTEXT_N_FFT       1024u
#define SPECTEXT_HOP           16u
    printf("Spectrogram text:\n");

    {
        const double st_dur = 2.25;
        const double st_pad = 0.5;   /* silence before and after text */
        unsigned st_text_max = (unsigned)(SPECTEXT_SAMPLE_RATE * st_dur) + 1024;
        double *st_text = malloc(st_text_max * sizeof(double));
        if (!st_text) {
            fprintf(stderr, "allocation failed for spectrogram text signal\n");
            free(conv_in);
            free(work_b);
            free(work_a);
            free(fx);
            free(buf);
            return 1;
        }
        unsigned st_text_len = MD_spectrogram_text(st_text, st_text_max, "HELLO",
                                                   400.0, 7300.0, st_dur,
                                                   SPECTEXT_SAMPLE_RATE);

        /* Pad with silence before and after */
        unsigned st_pad_samp = (unsigned)(st_pad * SPECTEXT_SAMPLE_RATE);
        unsigned st_len = st_pad_samp + st_text_len + st_pad_samp;
        double *st_sig = calloc(st_len, sizeof(double));
        if (!st_sig) {
            fprintf(stderr, "allocation failed for padded spectrogram text\n");
            free(st_text);
            free(conv_in);
            free(work_b);
            free(work_a);
            free(fx);
            free(buf);
            return 1;
        }
        memcpy(st_sig + st_pad_samp, st_text, st_text_len * sizeof(double));
        free(st_text);

        const unsigned st_bins   = SPECTEXT_N_FFT / 2 + 1;
        const unsigned st_frames = MD_stft_num_frames(st_len, SPECTEXT_N_FFT,
                                                      SPECTEXT_HOP);

        double *st_mag = malloc((size_t)st_frames * st_bins * sizeof(double));
        if (!st_mag) {
            fprintf(stderr, "allocation failed for spectrogram text STFT\n");
            free(st_sig);
            free(conv_in);
            free(work_b);
            free(work_a);
            free(fx);
            free(buf);
            return 1;
        }
        MD_stft(st_sig, st_len, SPECTEXT_N_FFT, SPECTEXT_HOP, st_mag);

        const char *st_path  = "guides/plots/spectext_hello_spectrogram.html";
        const char *st_title = "Spectrogram Text \"HELLO\"";
        FILE *fp = fopen(st_path, "w");
        if (!fp) {
            fprintf(stderr, "cannot open %s for writing\n", st_path);
            free(st_mag);
            free(st_sig);
            free(conv_in);
            free(work_b);
            free(work_a);
            free(fx);
            free(buf);
            return 1;
        }

        write_head(fp, st_title);

        /* Time axis */
        const double st_time_step = (double)SPECTEXT_HOP / SPECTEXT_SAMPLE_RATE;
        fprintf(fp,
                "    const times = Array.from({length: %u}, (_, f) => f * %.10g);\n",
                st_frames, st_time_step);

        /* Frequency axis */
        const double st_freq_step = SPECTEXT_SAMPLE_RATE / (double)SPECTEXT_N_FFT;
        fprintf(fp,
                "    const freqs = Array.from({length: %u}, (_, k) => k * %.10g);\n",
                st_bins, st_freq_step);

        /* Spectrogram matrix: z[k][f] — use "%.0f" for integer dB to reduce file size */
        fprintf(fp, "    const z = [\n");
        for (unsigned k = 0; k < st_bins; k++) {
            fprintf(fp, "      [");
            for (unsigned f = 0; f < st_frames; f++) {
                double db = 20.0 * log10(fmax(st_mag[f * st_bins + k]
                                              / (double)SPECTEXT_N_FFT, 1e-6));
                fprintf(fp, "%.0f", db);
                if (f + 1 < st_frames) fprintf(fp, ",");
            }
            fprintf(fp, "]");
            if (k + 1 < st_bins) fprintf(fp, ",");
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
            "      hovertemplate: 't: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>%%{z:.0f} dB<extra></extra>'\n"
            "    }], {\n"
            "      title: { text: '%s', font: { size: 13 } },\n"
            "      xaxis: { title: 'Time (s)' },\n"
            "      yaxis: { title: 'Frequency (Hz)', range: [0, 8000] },\n"
            "      margin: { t: 35, r: 80, b: 50, l: 60 },\n"
            "      height: 360\n"
            "    }, { responsive: true });\n",
            st_title);

        write_foot(fp);
        fclose(fp);
        printf("  %s\n", st_path);

        free(st_mag);
        free(st_sig);
    }

    /* ================================================================
     * Shepard tone spectrograms (rising + falling)
     * ================================================================*/
    {
        const double shep_sr = 44100.0;
        const double shep_dur = 5.0;
        const unsigned shep_n = (unsigned)(shep_sr * shep_dur);
        const unsigned shep_fft = 2048;
        const unsigned shep_hop = 512;

        double *shep_sig = malloc(shep_n * sizeof(double));
        if (!shep_sig) {
            fprintf(stderr, "allocation failed for shepard signal\n");
            free(conv_in); free(work_b); free(work_a); free(fx); free(buf);
            return 1;
        }

        struct { const char *path; const char *title; double rate; } shep_cfg[] = {
            { "guides/plots/shepard_rising_spectrogram.html",
              "Shepard Tone (Rising 0.5 oct/s)",  0.5 },
            { "guides/plots/shepard_falling_spectrogram.html",
              "Shepard Tone (Falling 0.5 oct/s)", -0.5 },
        };

        for (int sc = 0; sc < 2; sc++) {
            MD_shepard_tone(shep_sig, shep_n, 0.8, 440.0, shep_sr,
                            shep_cfg[sc].rate, 8);

            unsigned shep_frames = MD_stft_num_frames(shep_n, shep_fft, shep_hop);
            unsigned shep_bins   = shep_fft / 2 + 1;
            double *shep_mag = malloc((size_t)shep_frames * shep_bins * sizeof(double));
            if (!shep_mag) {
                fprintf(stderr, "allocation failed for shepard STFT\n");
                continue;
            }
            MD_stft(shep_sig, shep_n, shep_fft, shep_hop, shep_mag);

            FILE *fp = fopen(shep_cfg[sc].path, "w");
            if (!fp) {
                fprintf(stderr, "cannot open %s\n", shep_cfg[sc].path);
                free(shep_mag);
                continue;
            }

            write_head(fp, shep_cfg[sc].title);

            /* Time axis */
            fprintf(fp, "    const times = [");
            for (unsigned f = 0; f < shep_frames; f++) {
                fprintf(fp, "%.4f", (double)(f * shep_hop) / shep_sr);
                if (f + 1 < shep_frames) fprintf(fp, ",");
            }
            fprintf(fp, "];\n");

            /* Frequency axis */
            fprintf(fp, "    const freqs = [");
            for (unsigned k = 0; k < shep_bins; k++) {
                fprintf(fp, "%.2f", (double)k * shep_sr / (double)shep_fft);
                if (k + 1 < shep_bins) fprintf(fp, ",");
            }
            fprintf(fp, "];\n");

            /* Z matrix (dB) */
            fprintf(fp, "    const z = [\n");
            for (unsigned k = 0; k < shep_bins; k++) {
                fprintf(fp, "      [");
                for (unsigned f = 0; f < shep_frames; f++) {
                    double db = 20.0 * log10(fmax(shep_mag[f * shep_bins + k]
                                                   / (double)shep_fft, 1e-6));
                    fprintf(fp, "%.1f", db);
                    if (f + 1 < shep_frames) fprintf(fp, ",");
                }
                fprintf(fp, "]");
                if (k + 1 < shep_bins) fprintf(fp, ",");
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
                "      hovertemplate: 't: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>%%{z:.0f} dB<extra></extra>'\n"
                "    }], {\n"
                "      title: { text: '%s', font: { size: 13 } },\n"
                "      xaxis: { title: 'Time (s)' },\n"
                "      yaxis: { title: 'Frequency (Hz)', type: 'log', range: [%s, %s] },\n"
                "      margin: { t: 35, r: 80, b: 50, l: 60 },\n"
                "      height: 360\n"
                "    }, { responsive: true });\n",
                shep_cfg[sc].title,
                "Math.log10(30)", "Math.log10(10000)");

            write_foot(fp);
            fclose(fp);
            printf("  %s\n", shep_cfg[sc].path);

            free(shep_mag);
        }
        free(shep_sig);
    }

    /* ===================================================================
     * Steganography plots
     * =================================================================== */
    printf("\nSteganography plots:\n");
    {
        const double steg_sr = 44100.0;
        const unsigned steg_n = (unsigned)(steg_sr * 3.0);
        double *steg_host  = malloc(steg_n * sizeof(double));
        double *steg_lsb   = malloc(steg_n * sizeof(double));
        double *steg_freq  = malloc(steg_n * sizeof(double));
        if (!steg_host || !steg_lsb || !steg_freq) {
            fprintf(stderr, "allocation failed for steganography plots\n");
        } else {
            const char *secret = "Hidden message inside audio!";
            MD_sine_wave(steg_host, steg_n, 0.8, 440.0, steg_sr);
            MD_steg_encode(steg_host, steg_lsb, steg_n, steg_sr,
                           secret, MD_STEG_LSB);
            MD_steg_encode(steg_host, steg_freq, steg_n, steg_sr,
                           secret, MD_STEG_FREQ_BAND);

            /* --- LSB difference plot --- */
            {
                FILE *fp = fopen("guides/plots/steg_lsb_diff.html", "w");
                if (!fp) { fprintf(stderr, "cannot open steg_lsb_diff.html\n"); }
                else {
                    write_head(fp, "LSB Steganography — Difference Signal");

                    /* Compute difference (host - stego) for first 2000 samples */
                    unsigned show_n = 2000;
                    fprintf(fp, "    const t = [");
                    for (unsigned i = 0; i < show_n; i++) {
                        fprintf(fp, "%.6f", (double)i / steg_sr);
                        if (i + 1 < show_n) fprintf(fp, ",");
                    }
                    fprintf(fp, "];\n");

                    fprintf(fp, "    const diff = [");
                    for (unsigned i = 0; i < show_n; i++) {
                        fprintf(fp, "%.8e", steg_host[i] - steg_lsb[i]);
                        if (i + 1 < show_n) fprintf(fp, ",");
                    }
                    fprintf(fp, "];\n");

                    fprintf(fp,
                        "    Plotly.newPlot('plot', [{\n"
                        "      x: t, y: diff,\n"
                        "      type: 'scatter', mode: 'lines',\n"
                        "      line: { width: 1, color: '#e74c3c' },\n"
                        "      name: 'host - stego'\n"
                        "    }], {\n"
                        "      title: { text: 'LSB Difference (host \\u2212 stego)', font: { size: 13 } },\n"
                        "      xaxis: { title: 'Time (s)' },\n"
                        "      yaxis: { title: 'Amplitude', tickformat: '.1e' },\n"
                        "      margin: { t: 35, r: 30, b: 50, l: 70 },\n"
                        "      height: 360\n"
                        "    }, { responsive: true });\n");

                    write_foot(fp);
                    fclose(fp);
                    printf("  guides/plots/steg_lsb_diff.html\n");
                }
            }

            /* --- Frequency-band spectrogram --- */
            {
                /* Use a small portion (0.5 s) to show the BFSK carriers */
                const unsigned spec_dur_samp = (unsigned)(0.5 * steg_sr);
                const unsigned spec_fft = 2048u;
                const unsigned spec_hop = 64u;
                const unsigned spec_bins = spec_fft / 2 + 1;
                unsigned spec_frames = (spec_dur_samp >= spec_fft)
                    ? (spec_dur_samp - spec_fft) / spec_hop + 1 : 0;

                double *win = malloc(spec_fft * sizeof(double));
                double *frm = malloc(spec_fft * sizeof(double));
                double *mag = malloc(spec_bins * sizeof(double));
                double *all_mag = nullptr;
                if (spec_frames > 0)
                    all_mag = malloc(spec_frames * spec_bins * sizeof(double));

                if (win && frm && mag && all_mag) {
                    MD_Gen_Hann_Win(win, spec_fft);

                    for (unsigned f = 0; f < spec_frames; f++) {
                        unsigned start = f * spec_hop;
                        for (unsigned i = 0; i < spec_fft; i++)
                            frm[i] = steg_freq[start + i] * win[i];
                        MD_magnitude_spectrum(frm, spec_fft, mag);
                        for (unsigned k = 0; k < spec_bins; k++)
                            all_mag[f * spec_bins + k] = mag[k];
                    }

                    FILE *fp = fopen("guides/plots/steg_freq_spectrogram.html", "w");
                    if (!fp) {
                        fprintf(stderr, "cannot open steg_freq_spectrogram.html\n");
                    } else {
                        write_head(fp, "Frequency-Band Steganography — Spectrogram");

                        /* Times array */
                        fprintf(fp, "    const times = [");
                        for (unsigned f = 0; f < spec_frames; f++) {
                            fprintf(fp, "%.4f", (double)(f * spec_hop) / steg_sr);
                            if (f + 1 < spec_frames) fprintf(fp, ",");
                        }
                        fprintf(fp, "];\n");

                        /* Frequencies array */
                        fprintf(fp, "    const freqs = [");
                        for (unsigned k = 0; k < spec_bins; k++) {
                            fprintf(fp, "%.1f", (double)k * steg_sr / (double)spec_fft);
                            if (k + 1 < spec_bins) fprintf(fp, ",");
                        }
                        fprintf(fp, "];\n");

                        /* Z matrix (dB) */
                        fprintf(fp, "    const z = [\n");
                        for (unsigned k = 0; k < spec_bins; k++) {
                            fprintf(fp, "      [");
                            for (unsigned f = 0; f < spec_frames; f++) {
                                double db = 20.0 * log10(fmax(
                                    all_mag[f * spec_bins + k]
                                    / (double)spec_fft, 1e-6));
                                fprintf(fp, "%.1f", db);
                                if (f + 1 < spec_frames) fprintf(fp, ",");
                            }
                            fprintf(fp, "]");
                            if (k + 1 < spec_bins) fprintf(fp, ",");
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
                            "      hovertemplate: "
                            "'t: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>%%{z:.0f} dB<extra></extra>'\n"
                            "    }], {\n"
                            "      title: { text: 'Frequency-Band Stego (BFSK at 18.5/19.5 kHz)',"
                            " font: { size: 13 } },\n"
                            "      xaxis: { title: 'Time (s)' },\n"
                            "      yaxis: { title: 'Frequency (Hz)' },\n"
                            "      margin: { t: 35, r: 80, b: 50, l: 60 },\n"
                            "      height: 360\n"
                            "    }, { responsive: true });\n");

                        write_foot(fp);
                        fclose(fp);
                        printf("  guides/plots/steg_freq_spectrogram.html\n");
                    }
                }
                free(all_mag);
                free(mag);
                free(frm);
                free(win);
            }
        }
        free(steg_freq);
        free(steg_lsb);
        free(steg_host);
    }

    MD_shutdown();
    free(conv_in);
    free(work_b);
    free(work_a);
    free(fx);
    free(buf);
    return 0;
}
