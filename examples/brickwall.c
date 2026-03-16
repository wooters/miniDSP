/**
 * @file brickwall.c
 * @brief Example: FFT-based brickwall lowpass filter — before/after spectrum.
 *
 * Demonstrates MD_lowpass_brickwall() by:
 *   1. Generating a mixed signal (400 Hz + 12000 Hz) at 48 kHz.
 *   2. Computing the magnitude spectrum before filtering.
 *   3. Applying the brickwall lowpass at 4000 Hz.
 *   4. Computing the magnitude spectrum after filtering.
 *   5. Writing a two-subplot Plotly HTML with before/after spectra.
 *
 * Build and run (from the repository root):
 *   make -C examples brickwall
 *   cd examples && ./brickwall
 *   open examples/brickwall.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

/* -----------------------------------------------------------------------
 * Compact iframe-friendly HTML (no <h1> or subtitle).
 * -----------------------------------------------------------------------*/
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
        "  <div id=\"before-plot\"></div>\n"
        "  <div id=\"after-plot\"></div>\n"
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

int main(void)
{
    /* ------------------------------------------------------------------
     * Signal parameters
     * ----------------------------------------------------------------*/
    const unsigned   N           = 8192;
    const double     sample_rate = 48000.0;
    const unsigned   num_bins    = N / 2 + 1;

    const double     freq_lo     =   400.0;   /* below cutoff  */
    const double     freq_hi     = 12000.0;   /* above cutoff  */
    const double     cutoff_hz   =  4000.0;

    /* ------------------------------------------------------------------
     * Generate the mixed signal
     * ----------------------------------------------------------------*/
    //! [generate-signal]
    double *signal = malloc(N * sizeof(double));
    double *tone   = malloc(N * sizeof(double));

    MD_sine_wave(signal, N, 1.0, freq_lo, sample_rate);
    MD_sine_wave(tone,   N, 0.8, freq_hi, sample_rate);
    for (unsigned i = 0; i < N; i++)
        signal[i] += tone[i];

    free(tone);
    //! [generate-signal]

    /* ------------------------------------------------------------------
     * Compute magnitude spectrum BEFORE filtering
     * ----------------------------------------------------------------*/
    double *mag_before = malloc(num_bins * sizeof(double));
    double *freqs      = malloc(num_bins * sizeof(double));

    MD_magnitude_spectrum(signal, N, mag_before);

    for (unsigned k = 0; k < num_bins; k++) {
        freqs[k] = (double)k * sample_rate / (double)N;
        mag_before[k] /= (double)N;
        if (k > 0 && k < N / 2)
            mag_before[k] *= 2.0;
    }

    /* ------------------------------------------------------------------
     * Apply the brickwall lowpass filter
     * ----------------------------------------------------------------*/
    //! [apply-brickwall]
    MD_lowpass_brickwall(signal, N, cutoff_hz, sample_rate);
    //! [apply-brickwall]

    /* ------------------------------------------------------------------
     * Compute magnitude spectrum AFTER filtering
     * ----------------------------------------------------------------*/
    double *mag_after = malloc(num_bins * sizeof(double));

    MD_magnitude_spectrum(signal, N, mag_after);

    for (unsigned k = 0; k < num_bins; k++) {
        mag_after[k] /= (double)N;
        if (k > 0 && k < N / 2)
            mag_after[k] *= 2.0;
    }

    /* ------------------------------------------------------------------
     * Write HTML visualisation
     * ----------------------------------------------------------------*/
    const char *html_file = "brickwall.html";
    FILE *fp = fopen(html_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", html_file);
        return 1;
    }

    write_head(fp, "Brickwall Lowpass — Before / After");

    /* Embed spectrum data as JS arrays */
    plot_html_js_array(fp, "freqs", freqs, num_bins, "%.2f");
    plot_html_js_array(fp, "magBefore", mag_before, num_bins, "%.8f");
    plot_html_js_array(fp, "magAfter", mag_after, num_bins, "%.8f");

    fprintf(fp,
        "\n"
        "    var dbBefore = magBefore.map(m => 20*Math.log10(Math.max(m,1e-6)));\n"
        "    var dbAfter  = magAfter.map(m => 20*Math.log10(Math.max(m,1e-6)));\n"
        "\n"
        "    Plotly.newPlot('before-plot', [{\n"
        "      x: freqs, y: dbBefore,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#2563eb', width: 1.2 },\n"
        "      hovertemplate: '%%{x:.0f} Hz<br>%%{y:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: 'Before (%.0f Hz + %.0f Hz)', font: { size: 14 } },\n"
        "      xaxis: { title: 'Frequency (Hz)' },\n"
        "      yaxis: { title: 'Magnitude (dB)', range: [-100, 5] },\n"
        "      shapes: [{ type: 'line', x0: %.0f, x1: %.0f,\n"
        "                  y0: -100, y1: 5,\n"
        "                  line: { color: '#dc2626', width: 1.5, dash: 'dash' } }],\n"
        "      annotations: [{ x: %.0f, y: 3, text: 'cutoff',\n"
        "                      showarrow: false, font: { size: 11, color: '#dc2626' } }],\n"
        "      margin: { t: 40, r: 20, b: 45, l: 55 },\n"
        "      height: 280\n"
        "    }, { responsive: true });\n"
        "\n"
        "    Plotly.newPlot('after-plot', [{\n"
        "      x: freqs, y: dbAfter,\n"
        "      type: 'scatter', mode: 'lines',\n"
        "      line: { color: '#16a34a', width: 1.2 },\n"
        "      hovertemplate: '%%{x:.0f} Hz<br>%%{y:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      title: { text: 'After — brickwall at %.0f Hz', font: { size: 14 } },\n"
        "      xaxis: { title: 'Frequency (Hz)' },\n"
        "      yaxis: { title: 'Magnitude (dB)', range: [-100, 5] },\n"
        "      shapes: [{ type: 'line', x0: %.0f, x1: %.0f,\n"
        "                  y0: -100, y1: 5,\n"
        "                  line: { color: '#dc2626', width: 1.5, dash: 'dash' } }],\n"
        "      annotations: [{ x: %.0f, y: 3, text: 'cutoff',\n"
        "                      showarrow: false, font: { size: 11, color: '#dc2626' } }],\n"
        "      margin: { t: 40, r: 20, b: 45, l: 55 },\n"
        "      height: 280\n"
        "    }, { responsive: true });\n",
        freq_lo, freq_hi,
        cutoff_hz, cutoff_hz, cutoff_hz,
        cutoff_hz,
        cutoff_hz, cutoff_hz, cutoff_hz);

    write_foot(fp);
    fclose(fp);

    printf("Wrote %s  (%.0f Hz + %.0f Hz, brickwall at %.0f Hz)\n",
           html_file, freq_lo, freq_hi, cutoff_hz);

    /* ------------------------------------------------------------------
     * Cleanup
     * ----------------------------------------------------------------*/
    free(mag_after);
    free(freqs);
    free(mag_before);
    free(signal);
    MD_shutdown();

    return 0;
}
