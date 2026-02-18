/**
 * @file basic_operations.c
 * @brief Example: basic signal analysis operations.
 *
 * This program demonstrates five fundamental signal analysis functions:
 *   1. MD_rms() -- root mean square (signal loudness)
 *   2. MD_zero_crossing_rate() -- sign-change frequency
 *   3. MD_autocorrelation() -- self-similarity at different lags
 *   4. MD_peak_detect() -- finding local maxima
 *   5. MD_mix() -- weighted sum of two signals
 *
 * It generates a sine wave and white noise, mixes them, then computes
 * measurements on each.  Results are written to CSV and an interactive
 * HTML visualisation (Plotly.js).
 *
 * Build and run (from the repository root):
 *   make -C examples basic_operations
 *   cd examples && ./basic_operations
 *   open basic_operations.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

int main(void)
{
    /* ------------------------------------------------------------------
     * Signal parameters
     * ----------------------------------------------------------------*/
    const unsigned N           = 4096;
    const double   sample_rate = 8000.0;
    const double   freq_hz     = 200.0;
    const double   amplitude   = 1.0;
    const unsigned max_lag     = 200;

    /* ------------------------------------------------------------------
     * Generate signals
     * ----------------------------------------------------------------*/
    double *sine  = malloc(N * sizeof(double));
    double *noise = malloc(N * sizeof(double));
    double *mixed = malloc(N * sizeof(double));

    if (!sine || !noise || !mixed) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    //! [generate-signals]
    MD_sine_wave(sine, N, amplitude, freq_hz, sample_rate);
    MD_white_noise(noise, N, 0.3, 42);
    //! [generate-signals]

    /* ------------------------------------------------------------------
     * Mix: 80% sine + 20% noise
     * ----------------------------------------------------------------*/
    //! [mix-signals]
    MD_mix(sine, noise, mixed, N, 0.8, 0.2);
    //! [mix-signals]

    /* ------------------------------------------------------------------
     * RMS measurements
     * ----------------------------------------------------------------*/
    //! [rms-measurements]
    double rms_sine  = MD_rms(sine, N);
    double rms_noise = MD_rms(noise, N);
    double rms_mixed = MD_rms(mixed, N);
    //! [rms-measurements]

    printf("RMS: sine=%.4f  noise=%.4f  mixed=%.4f\n",
           rms_sine, rms_noise, rms_mixed);

    /* ------------------------------------------------------------------
     * Zero-crossing rate
     * ----------------------------------------------------------------*/
    //! [zcr-measurements]
    double zcr_sine  = MD_zero_crossing_rate(sine, N);
    double zcr_noise = MD_zero_crossing_rate(noise, N);
    double zcr_mixed = MD_zero_crossing_rate(mixed, N);
    //! [zcr-measurements]

    printf("ZCR: sine=%.4f  noise=%.4f  mixed=%.4f\n",
           zcr_sine, zcr_noise, zcr_mixed);

    /* ------------------------------------------------------------------
     * Autocorrelation of the mixed signal
     * ----------------------------------------------------------------*/
    double *acf = malloc(max_lag * sizeof(double));
    if (!acf) { fprintf(stderr, "allocation failed\n"); return 1; }

    //! [autocorrelation]
    MD_autocorrelation(mixed, N, acf, max_lag);
    //! [autocorrelation]

    /* ------------------------------------------------------------------
     * Peak detection on the autocorrelation
     * ----------------------------------------------------------------*/
    unsigned *peaks = malloc(max_lag * sizeof(unsigned));
    unsigned num_peaks = 0;
    if (!peaks) { fprintf(stderr, "allocation failed\n"); return 1; }

    //! [peak-detection]
    MD_peak_detect(acf, max_lag, 0.3, 10, peaks, &num_peaks);
    //! [peak-detection]

    printf("Autocorrelation peaks (lag, value):");
    for (unsigned i = 0; i < num_peaks; i++) {
        printf(" (%u, %.3f)", peaks[i], acf[peaks[i]]);
    }
    printf("\n");

    if (num_peaks > 0) {
        double estimated_freq = sample_rate / (double)peaks[0];
        printf("Estimated fundamental frequency: %.1f Hz (true: %.1f Hz)\n",
               estimated_freq, freq_hz);
    }

    /* ------------------------------------------------------------------
     * Write CSV
     * ----------------------------------------------------------------*/
    const char *csv_file = "basic_operations.csv";
    FILE *fp = fopen(csv_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", csv_file);
        return 1;
    }

    /* First section: time-domain signals (first 400 samples) */
    unsigned show_n = 400;
    fprintf(fp, "sample,sine,noise,mixed\n");
    for (unsigned i = 0; i < show_n && i < N; i++) {
        fprintf(fp, "%u,%.6f,%.6f,%.6f\n", i, sine[i], noise[i], mixed[i]);
    }

    fclose(fp);
    printf("Wrote %u samples to %s\n", show_n, csv_file);

    /* ------------------------------------------------------------------
     * Write interactive HTML visualisation
     * ----------------------------------------------------------------*/
    const char *html_file = "basic_operations.html";
    fp = fopen(html_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", html_file);
        return 1;
    }

    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
             "%.0f Hz sine + noise &nbsp;|&nbsp; sample rate %.0f Hz "
             "&nbsp;|&nbsp; N = %u",
             freq_hz, sample_rate, N);

    plot_html_begin(fp, "Basic Signal Operations", subtitle, 0);

    /* Three plot containers */
    fprintf(fp,
        "  <div id=\"signals\" class=\"plot-container\"></div>\n"
        "  <div id=\"acf\" class=\"plot-container\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>What you are seeing:</strong><br>\n"
        "    <b>Top:</b> A pure sine wave, white noise, and their weighted mix "
        "(0.8&times;sine + 0.2&times;noise).<br>\n"
        "    <b>Bottom:</b> The autocorrelation of the mixed signal. "
        "Peaks at multiples of the fundamental period (%.1f samples = 1/%.0f Hz) "
        "reveal the dominant frequency despite the noise.<br>\n"
        "    <br>\n"
        "    <b>RMS:</b> sine=%.4f &nbsp; noise=%.4f &nbsp; mixed=%.4f<br>\n"
        "    <b>ZCR:</b> sine=%.4f &nbsp; noise=%.4f &nbsp; mixed=%.4f\n"
        "  </div>\n\n",
        sample_rate / freq_hz, freq_hz,
        rms_sine, rms_noise, rms_mixed,
        zcr_sine, zcr_noise, zcr_mixed);

    /* Embed data as JS arrays */
    fprintf(fp, "  <script>\n");

    /* Time-domain signals (first 400 samples) */
    double *samples_idx = malloc(show_n * sizeof(double));
    for (unsigned i = 0; i < show_n; i++) samples_idx[i] = (double)i;
    plot_html_js_array(fp, "idx", samples_idx, show_n, "%.0f");

    /* Re-use the signal arrays but only the first show_n elements */
    plot_html_js_array(fp, "sineData", sine, show_n, "%.6f");
    plot_html_js_array(fp, "noiseData", noise, show_n, "%.6f");
    plot_html_js_array(fp, "mixedData", mixed, show_n, "%.6f");

    /* Autocorrelation */
    double *lag_idx = malloc(max_lag * sizeof(double));
    for (unsigned i = 0; i < max_lag; i++) lag_idx[i] = (double)i;
    plot_html_js_array(fp, "lags", lag_idx, max_lag, "%.0f");
    plot_html_js_array(fp, "acfData", acf, max_lag, "%.6f");

    /* Peak markers */
    fprintf(fp, "    const peakLags = [");
    for (unsigned i = 0; i < num_peaks; i++) {
        fprintf(fp, "%u%s", peaks[i], i + 1 < num_peaks ? "," : "");
    }
    fprintf(fp, "];\n");
    fprintf(fp, "    const peakVals = [");
    for (unsigned i = 0; i < num_peaks; i++) {
        fprintf(fp, "%.6f%s", acf[peaks[i]], i + 1 < num_peaks ? "," : "");
    }
    fprintf(fp, "];\n\n");

    /* Plotly: signals subplot */
    fprintf(fp,
        "    Plotly.newPlot('signals', [\n"
        "      { x: idx, y: sineData, name: 'Sine', type: 'scatter',\n"
        "        line: { color: '#2563eb', width: 1 } },\n"
        "      { x: idx, y: noiseData, name: 'Noise', type: 'scatter',\n"
        "        line: { color: '#dc2626', width: 1 } },\n"
        "      { x: idx, y: mixedData, name: 'Mixed', type: 'scatter',\n"
        "        line: { color: '#059669', width: 1.2 } }\n"
        "    ], {\n"
        "      title: { text: 'Time-Domain Signals (first %u samples)', font: { size: 15 } },\n"
        "      xaxis: { title: 'Sample' },\n"
        "      yaxis: { title: 'Amplitude' },\n"
        "      margin: { t: 55, r: 30, b: 55, l: 65 },\n"
        "      height: 350,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n\n",
        show_n);

    /* Plotly: autocorrelation plot with peak markers */
    fprintf(fp,
        "    Plotly.newPlot('acf', [\n"
        "      { x: lags, y: acfData, name: 'Autocorrelation', type: 'scatter',\n"
        "        line: { color: '#7c3aed', width: 1.2 } },\n"
        "      { x: peakLags, y: peakVals, name: 'Peaks', mode: 'markers',\n"
        "        marker: { color: '#dc2626', size: 8, symbol: 'diamond' } }\n"
        "    ], {\n"
        "      title: { text: 'Autocorrelation of Mixed Signal', font: { size: 15 } },\n"
        "      xaxis: { title: 'Lag (samples)' },\n"
        "      yaxis: { title: 'Normalised R[lag]', range: [-1.1, 1.1] },\n"
        "      margin: { t: 55, r: 30, b: 55, l: 65 },\n"
        "      height: 350,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n");

    fprintf(fp, "  </script>\n");
    plot_html_end(fp);
    fclose(fp);
    printf("Wrote interactive plot to %s\n", html_file);

    /* ------------------------------------------------------------------
     * Cleanup (no MD_shutdown needed -- these functions are stateless)
     * ----------------------------------------------------------------*/
    free(samples_idx);
    free(lag_idx);
    free(peaks);
    free(acf);
    free(mixed);
    free(noise);
    free(sine);

    return 0;
}
