/**
 * @file pitch_detection.c
 * @brief Example: frame-based F0 tracking via autocorrelation and FFT peak-picking.
 *
 * This example generates a synthetic voiced signal with changing pitch,
 * then estimates F0 per frame using:
 *   1) MD_f0_autocorrelation()
 *   2) MD_f0_fft()
 *
 * Results are written to:
 *   - pitch_detection.csv
 *   - pitch_detection.html (interactive Plotly visualisation)
 *
 * Build and run (from the repository root):
 *   make -C examples pitch_detection
 *   cd examples && ./pitch_detection
 *   open pitch_detection.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

static int write_html(const char *path,
                      const double *times,
                      const double *f0_true,
                      const double *f0_acf,
                      const double *f0_fft,
                      unsigned num_frames,
                      double fs, unsigned frame_len, unsigned hop)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
             "Frame-based F0 tracking &nbsp;|&nbsp; fs=%.0f Hz "
             "&nbsp;|&nbsp; frame=%u &nbsp;|&nbsp; hop=%u",
             fs, frame_len, hop);

    plot_html_begin(fp, "Pitch Detection (F0 Tracking)", subtitle, 0);

    fprintf(fp,
        "  <div id=\"f0-track\" class=\"plot-container\"></div>\n"
        "  <div id=\"f0-error\" class=\"plot-container\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>What you are seeing:</strong><br>\n"
        "    <b>Top:</b> Ground-truth F0 and frame-wise estimates from "
        "autocorrelation and FFT peak-picking.<br>\n"
        "    <b>Bottom:</b> Absolute estimation error (Hz) per frame.<br>\n"
        "    ACF usually tracks fundamentals more reliably when harmonics/noise are present.\n"
        "  </div>\n");

    fprintf(fp, "  <script>\n");
    plot_html_js_array(fp, "times", times, num_frames, "%.6f");
    plot_html_js_array(fp, "f0True", f0_true, num_frames, "%.6f");
    plot_html_js_array(fp, "f0Acf", f0_acf, num_frames, "%.6f");
    plot_html_js_array(fp, "f0Fft", f0_fft, num_frames, "%.6f");

    fprintf(fp,
        "    const errAcf = f0Acf.map((v, i) => (v > 0 ? Math.abs(v - f0True[i]) : null));\n"
        "    const errFft = f0Fft.map((v, i) => (v > 0 ? Math.abs(v - f0True[i]) : null));\n"
        "\n"
        "    Plotly.newPlot('f0-track', [\n"
        "      { x: times, y: f0True, name: 'Ground Truth', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#111827', width: 2 } },\n"
        "      { x: times, y: f0Acf, name: 'Autocorrelation F0', type: 'scatter', mode: 'lines+markers',\n"
        "        marker: { size: 4 }, line: { color: '#2563eb', width: 1.4 } },\n"
        "      { x: times, y: f0Fft, name: 'FFT F0', type: 'scatter', mode: 'lines+markers',\n"
        "        marker: { size: 4 }, line: { color: '#dc2626', width: 1.2 } }\n"
        "    ], {\n"
        "      title: { text: 'F0 Tracks vs Ground Truth', font: { size: 15 } },\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'Frequency (Hz)' },\n"
        "      margin: { t: 55, r: 30, b: 55, l: 65 },\n"
        "      height: 360,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n"
        "\n"
        "    Plotly.newPlot('f0-error', [\n"
        "      { x: times, y: errAcf, name: '|ACF - truth|', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#2563eb', width: 1.4 } },\n"
        "      { x: times, y: errFft, name: '|FFT - truth|', type: 'scatter', mode: 'lines',\n"
        "        line: { color: '#dc2626', width: 1.3 } }\n"
        "    ], {\n"
        "      title: { text: 'Absolute F0 Error per Frame', font: { size: 15 } },\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'Absolute Error (Hz)', rangemode: 'tozero' },\n"
        "      margin: { t: 55, r: 30, b: 55, l: 65 },\n"
        "      height: 320,\n"
        "      legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "    }, { responsive: true });\n"
        "  </script>\n");

    plot_html_end(fp);
    fclose(fp);
    return 0;
}

static double mean_abs_error(const double *est, const double *truth, unsigned N)
{
    double sum = 0.0;
    unsigned count = 0;
    for (unsigned i = 0; i < N; i++) {
        if (est[i] > 0.0) {
            sum += fabs(est[i] - truth[i]);
            count++;
        }
    }
    return (count > 0) ? (sum / (double)count) : 0.0;
}

int main(void)
{
    const double sample_rate = 16000.0;
    const double duration_s  = 2.0;
    const unsigned N = (unsigned)(sample_rate * duration_s);

    const unsigned frame_len = 1024;
    const unsigned hop = 256;
    const double min_f0 = 80.0;
    const double max_f0 = 400.0;

    const unsigned num_frames = (N - frame_len) / hop + 1;

    double *signal   = malloc(N * sizeof(double));
    double *noise    = malloc(N * sizeof(double));
    double *f0_true  = malloc(num_frames * sizeof(double));
    double *f0_acf   = malloc(num_frames * sizeof(double));
    double *f0_fft   = malloc(num_frames * sizeof(double));
    double *times    = malloc(num_frames * sizeof(double));

    if (!signal || !noise || !f0_true || !f0_acf || !f0_fft || !times) {
        fprintf(stderr, "allocation failed\n");
        free(times);
        free(f0_fft);
        free(f0_acf);
        free(f0_true);
        free(noise);
        free(signal);
        return 1;
    }

    /* ------------------------------------------------------------------
     * Build a synthetic voiced signal with piecewise-constant F0.
     * ----------------------------------------------------------------*/
    const unsigned seg1 = N / 3;
    const unsigned seg2 = 2 * N / 3;
    double phase = 0.0;
    MD_white_noise(noise, N, 0.08, 42);

    for (unsigned n = 0; n < N; n++) {
        double f0 = (n < seg1) ? 140.0 : (n < seg2 ? 220.0 : 320.0);
        double dphi = 2.0 * M_PI * f0 / sample_rate;
        phase += dphi;

        /* Fundamental + harmonics + light noise */
        signal[n] = 0.75 * sin(phase)
                  + 0.22 * sin(2.0 * phase)
                  + 0.12 * sin(3.0 * phase)
                  + noise[n];
    }

    /* ------------------------------------------------------------------
     * Frame-wise F0 tracking
     * ----------------------------------------------------------------*/
    //! [frame-tracking]
    for (unsigned f = 0; f < num_frames; f++) {
        unsigned start = f * hop;
        const double *frame = signal + start;
        unsigned center = start + frame_len / 2;
        if (center >= N) center = N - 1;

        /* Ground truth for this frame (piecewise-constant by construction). */
        f0_true[f] = (center < seg1) ? 140.0 : (center < seg2 ? 220.0 : 320.0);
        times[f] = (double)center / sample_rate;

        //! [acf-reading-algorithm]
        f0_acf[f] = MD_f0_autocorrelation(frame, frame_len, sample_rate,
                                          min_f0, max_f0);
        //! [acf-reading-algorithm]

        //! [fft-reading-algorithm]
        f0_fft[f] = MD_f0_fft(frame, frame_len, sample_rate,
                              min_f0, max_f0);
        //! [fft-reading-algorithm]
    }
    //! [frame-tracking]

    double mae_acf = mean_abs_error(f0_acf, f0_true, num_frames);
    double mae_fft = mean_abs_error(f0_fft, f0_true, num_frames);

    printf("Frames: %u  (frame_len=%u, hop=%u)\n", num_frames, frame_len, hop);
    printf("Search range: %.1f-%.1f Hz\n", min_f0, max_f0);
    printf("MAE: autocorrelation=%.2f Hz, fft=%.2f Hz\n", mae_acf, mae_fft);

    const char *csv_file = "pitch_detection.csv";
    FILE *fp = fopen(csv_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", csv_file);
        free(times);
        free(f0_fft);
        free(f0_acf);
        free(f0_true);
        free(noise);
        free(signal);
        return 1;
    }

    fprintf(fp, "frame,time_s,f0_true_hz,f0_acf_hz,f0_fft_hz\n");
    for (unsigned f = 0; f < num_frames; f++) {
        fprintf(fp, "%u,%.6f,%.6f,%.6f,%.6f\n",
                f, times[f], f0_true[f], f0_acf[f], f0_fft[f]);
    }
    fclose(fp);
    printf("Wrote frame-level tracks to %s\n", csv_file);

    const char *html_file = "pitch_detection.html";
    if (write_html(html_file, times, f0_true, f0_acf, f0_fft,
                   num_frames, sample_rate, frame_len, hop) == 0) {
        printf("Wrote interactive plot to %s\n", html_file);
    }

    MD_shutdown();
    free(times);
    free(f0_fft);
    free(f0_acf);
    free(f0_true);
    free(noise);
    free(signal);
    return 0;
}
