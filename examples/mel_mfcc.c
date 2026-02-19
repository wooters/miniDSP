/**
 * @file mel_mfcc.c
 * @brief Example: compute mel energies and MFCCs from a single frame.
 *
 * Outputs:
 *   - mel_mfcc.csv
 *   - mel_mfcc.html (interactive Plotly visualisation)
 *
 * Build and run (from repository root):
 *   make -C examples mel_mfcc
 *   cd examples && ./mel_mfcc
 *   open mel_mfcc.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

static int write_html(const char *path,
                      const double *mel_centers,
                      const double *mel_energies,
                      unsigned num_mels,
                      const double *mfcc,
                      unsigned num_coeffs,
                      unsigned N, double sample_rate)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
             "Single-frame analysis &nbsp;|&nbsp; N=%u &nbsp;|&nbsp; fs=%.0f Hz",
             N, sample_rate);
    plot_html_begin(fp, "Mel Filterbank and MFCC", subtitle, 0);

    fprintf(fp,
        "  <div id=\"mel-plot\" class=\"plot-container\"></div>\n"
        "  <div id=\"mfcc-plot\" class=\"plot-container\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>Top:</strong> mel-band energies (power-domain integration over triangular filters).<br>\n"
        "    <strong>Bottom:</strong> MFCC vector (DCT-II of log mel energies, C0 included).\n"
        "  </div>\n");

    fprintf(fp, "  <script>\n");
    plot_html_js_array(fp, "melCenters", mel_centers, num_mels, "%.6f");
    plot_html_js_array(fp, "melEnergies", mel_energies, num_mels, "%.10f");
    plot_html_js_array(fp, "mfccVals", mfcc, num_coeffs, "%.10f");
    fprintf(fp,
        "    const melIdx = Array.from({length: %u}, (_, i) => i);\n"
        "    const mfccIdx = Array.from({length: %u}, (_, i) => i);\n"
        "\n"
        "    Plotly.newPlot('mel-plot', [{\n"
        "      x: melCenters, y: melEnergies,\n"
        "      type: 'bar', marker: { color: '#2563eb' },\n"
        "      hovertemplate: 'Center: %%{x:.1f} Hz<br>Energy: %%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: 'Mel-Band Energies', font: { size: 15 } },\n"
        "      xaxis: { title: 'Approximate Band Center (Hz)' },\n"
        "      yaxis: { title: 'Energy', rangemode: 'tozero' },\n"
        "      margin: { t: 50, r: 30, b: 55, l: 70 },\n"
        "      height: 320\n"
        "    }, { responsive: true });\n"
        "\n"
        "    Plotly.newPlot('mfcc-plot', [{\n"
        "      x: mfccIdx, y: mfccVals,\n"
        "      type: 'bar', marker: { color: '#dc2626' },\n"
        "      hovertemplate: 'c%%{x}: %%{y:.6f}<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: 'MFCCs (C0 included)', font: { size: 15 } },\n"
        "      xaxis: { title: 'Coefficient Index' },\n"
        "      yaxis: { title: 'Value' },\n"
        "      margin: { t: 50, r: 30, b: 55, l: 70 },\n"
        "      height: 320\n"
        "    }, { responsive: true });\n"
        "  </script>\n",
        num_mels, num_coeffs);

    plot_html_end(fp);
    fclose(fp);
    return 0;
}

int main(void)
{
    const unsigned N = 1024;
    const double sample_rate = 16000.0;
    const unsigned num_bins = N / 2 + 1;
    const unsigned num_mels = 26;
    const unsigned num_coeffs = 13;
    const double min_freq_hz = 80.0;
    const double max_freq_hz = 7600.0;

    double *signal = malloc(N * sizeof(double));
    double *mel = malloc(num_mels * sizeof(double));
    double *mfcc = malloc(num_coeffs * sizeof(double));
    double *fb = malloc((size_t)num_mels * num_bins * sizeof(double));
    double *mel_centers = malloc(num_mels * sizeof(double));
    if (!signal || !mel || !mfcc || !fb || !mel_centers) {
        fprintf(stderr, "allocation failed\n");
        free(mel_centers);
        free(fb);
        free(mfcc);
        free(mel);
        free(signal);
        return 1;
    }

    for (unsigned n = 0; n < N; n++) {
        double t = (double)n / sample_rate;
        signal[n] = 0.7 * sin(2.0 * M_PI * 440.0 * t)
                  + 0.2 * cos(2.0 * M_PI * 1000.0 * t)
                  + 0.1 * sin(2.0 * M_PI * 3000.0 * t);
    }

    //! [compute-mel]
    MD_mel_energies(signal, N, sample_rate, num_mels,
                    min_freq_hz, max_freq_hz, mel);
    //! [compute-mel]

    //! [compute-mfcc]
    MD_mfcc(signal, N, sample_rate, num_mels, num_coeffs,
            min_freq_hz, max_freq_hz, mfcc);
    //! [compute-mfcc]

    MD_mel_filterbank(N, sample_rate, num_mels, min_freq_hz, max_freq_hz, fb);
    for (unsigned m = 0; m < num_mels; m++) {
        const double *row = fb + (size_t)m * num_bins;
        double sum_w = 0.0;
        double sum_fw = 0.0;
        for (unsigned k = 0; k < num_bins; k++) {
            double f_hz = (double)k * sample_rate / (double)N;
            sum_w += row[k];
            sum_fw += f_hz * row[k];
        }
        mel_centers[m] = (sum_w > 0.0) ? (sum_fw / sum_w) : 0.0;
    }

    FILE *csv = fopen("mel_mfcc.csv", "w");
    if (!csv) {
        fprintf(stderr, "cannot open mel_mfcc.csv for writing\n");
        free(mel_centers);
        free(fb);
        free(mfcc);
        free(mel);
        free(signal);
        return 1;
    }

    fprintf(csv, "kind,index,frequency_hz,value\n");
    for (unsigned m = 0; m < num_mels; m++) {
        fprintf(csv, "mel,%u,%.6f,%.10f\n", m, mel_centers[m], mel[m]);
    }
    for (unsigned c = 0; c < num_coeffs; c++) {
        fprintf(csv, "mfcc,%u,0.0,%.10f\n", c, mfcc[c]);
    }
    fclose(csv);
    printf("Wrote mel energies and MFCCs to mel_mfcc.csv\n");

    if (write_html("mel_mfcc.html", mel_centers, mel, num_mels,
                   mfcc, num_coeffs, N, sample_rate) == 0) {
        printf("Wrote interactive plot to mel_mfcc.html\n");
    }

    MD_shutdown();
    free(mel_centers);
    free(fb);
    free(mfcc);
    free(mel);
    free(signal);
    return 0;
}

