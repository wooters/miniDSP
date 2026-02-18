/**
 * @file simple_effects.c
 * @brief Example: simple effects (delay/echo, tremolo, comb reverb).
 *
 * Build and run (from repo root):
 *   make -C examples simple_effects
 *   ./examples/simple_effects
 *   open examples/simple_effects.html
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

int main(void)
{
    const unsigned N = 4096;
    const double sample_rate = 8000.0;

    double *delay_src = calloc(N, sizeof(double));
    double *delay_out = calloc(N, sizeof(double));
    double *trem_src  = calloc(N, sizeof(double));
    double *trem_out  = calloc(N, sizeof(double));
    double *comb_src  = calloc(N, sizeof(double));
    double *comb_out  = calloc(N, sizeof(double));

    if (!delay_src || !delay_out || !trem_src || !trem_out || !comb_src || !comb_out) {
        fprintf(stderr, "allocation failed\n");
        free(comb_out); free(comb_src); free(trem_out); free(trem_src); free(delay_out); free(delay_src);
        return 1;
    }

    /* Percussive click train source for delay/echo. */
    for (unsigned i = 0; i < N; i += 600) {
        delay_src[i] = 1.0;
    }

    /* Steady tone source for tremolo. */
    MD_sine_wave(trem_src, N, 0.8, 220.0, sample_rate);

    /* Decaying burst source for comb reverb. */
    for (unsigned i = 0; i < N; i++) {
        double env = exp(-(double)i / 800.0);
        comb_src[i] = 0.8 * env * sin(2.0 * M_PI * 330.0 * (double)i / sample_rate);
    }

    //! [delay-echo]
    MD_delay_echo(delay_src, delay_out, N, 400, 0.45, 1.0, 0.6);
    //! [delay-echo]

    //! [tremolo]
    MD_tremolo(trem_src, trem_out, N, 5.0, 0.8, sample_rate);
    //! [tremolo]

    //! [comb-reverb]
    MD_comb_reverb(comb_src, comb_out, N, 120, 0.75, 0.7, 0.6);
    //! [comb-reverb]

    printf("Delay RMS:  before=%.4f after=%.4f\n", MD_rms(delay_src, N), MD_rms(delay_out, N));
    printf("Tremolo RMS: before=%.4f after=%.4f\n", MD_rms(trem_src, N), MD_rms(trem_out, N));
    printf("Comb RMS:   before=%.4f after=%.4f\n", MD_rms(comb_src, N), MD_rms(comb_out, N));

    const unsigned show_n = 1200;
    FILE *fp = fopen("simple_effects.csv", "w");
    if (!fp) {
        fprintf(stderr, "cannot open simple_effects.csv for writing\n");
        free(comb_out); free(comb_src); free(trem_out); free(trem_src); free(delay_out); free(delay_src);
        return 1;
    }

    fprintf(fp, "n,delay_before,delay_after,trem_before,trem_after,comb_before,comb_after\n");
    for (unsigned i = 0; i < show_n && i < N; i++) {
        fprintf(fp, "%u,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                i, delay_src[i], delay_out[i], trem_src[i], trem_out[i], comb_src[i], comb_out[i]);
    }
    fclose(fp);

    fp = fopen("simple_effects.html", "w");
    if (!fp) {
        fprintf(stderr, "cannot open simple_effects.html for writing\n");
        free(comb_out); free(comb_src); free(trem_out); free(trem_src); free(delay_out); free(delay_src);
        return 1;
    }

    plot_html_begin(fp, "Simple Effects",
                    "Delay/Echo, Tremolo, and Comb Reverb (before vs after)", 0);

    fprintf(fp,
        "  <div id=\"delayPlot\" class=\"plot-container\"></div>\n"
        "  <div id=\"tremPlot\" class=\"plot-container\"></div>\n"
        "  <div id=\"combPlot\" class=\"plot-container\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>What you are seeing:</strong><br>\n"
        "    Delay: click train with repeating echoes.<br>\n"
        "    Tremolo: low-frequency amplitude modulation of a steady tone.<br>\n"
        "    Comb reverb: resonant, decaying echoes caused by feedback delay.\n"
        "  </div>\n");

    fprintf(fp, "  <script>\n");
    double *idx = malloc(show_n * sizeof(double));
    if (!idx) {
        fprintf(stderr, "allocation failed\n");
        fclose(fp);
        free(comb_out); free(comb_src); free(trem_out); free(trem_src); free(delay_out); free(delay_src);
        return 1;
    }
    for (unsigned i = 0; i < show_n; i++) idx[i] = (double)i;

    plot_html_js_array(fp, "idx", idx, show_n, "%.0f");
    plot_html_js_array(fp, "delayBefore", delay_src, show_n, "%.6f");
    plot_html_js_array(fp, "delayAfter", delay_out, show_n, "%.6f");
    plot_html_js_array(fp, "tremBefore", trem_src, show_n, "%.6f");
    plot_html_js_array(fp, "tremAfter", trem_out, show_n, "%.6f");
    plot_html_js_array(fp, "combBefore", comb_src, show_n, "%.6f");
    plot_html_js_array(fp, "combAfter", comb_out, show_n, "%.6f");

    fprintf(fp,
        "    function twoTrace(div, title, a, b) {\n"
        "      Plotly.newPlot(div, [\n"
        "        { x: idx, y: a, name: 'Before', type: 'scatter', line: { color: '#2563eb', width: 1.1 } },\n"
        "        { x: idx, y: b, name: 'After',  type: 'scatter', line: { color: '#dc2626', width: 1.1 } }\n"
        "      ], {\n"
        "        title: { text: title, font: { size: 15 } },\n"
        "        xaxis: { title: 'Sample Index' },\n"
        "        yaxis: { title: 'Amplitude' },\n"
        "        margin: { t: 50, r: 25, b: 50, l: 60 },\n"
        "        height: 300,\n"
        "        legend: { x: 1, y: 1, xanchor: 'right' }\n"
        "      }, { responsive: true });\n"
        "    }\n"
        "    twoTrace('delayPlot', 'Delay / Echo', delayBefore, delayAfter);\n"
        "    twoTrace('tremPlot', 'Tremolo', tremBefore, tremAfter);\n"
        "    twoTrace('combPlot', 'Comb Reverb', combBefore, combAfter);\n");

    fprintf(fp, "  </script>\n");
    plot_html_end(fp);
    fclose(fp);

    printf("Wrote simple_effects.csv and simple_effects.html\n");

    free(idx);
    free(comb_out);
    free(comb_src);
    free(trem_out);
    free(trem_src);
    free(delay_out);
    free(delay_src);
    return 0;
}
