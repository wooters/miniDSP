/**
 * @file spectrogram.c
 * @brief Example: compute and visualise the STFT spectrogram of a chirp signal.
 *
 * This program demonstrates MD_stft() by:
 *   1. Generating a linear chirp sweeping from 200 Hz to 4000 Hz over 2 s.
 *   2. Computing the STFT with a 32 ms Hanning window and 75% overlap.
 *   3. Converting the magnitudes to a dB scale.
 *   4. Writing results to CSV and to an interactive HTML heatmap (Plotly.js).
 *
 * A chirp (swept sine) is ideal for testing an STFT because its instantaneous
 * frequency changes continuously, producing a clearly visible diagonal stripe
 * across the time-frequency plane.
 *
 * Build and run (from the repository root):
 *   make                          # build libminidsp.a first
 *   make -C examples plot         # build example, run it, generate HTML
 *   open examples/spectrogram.html   # interactive spectrogram heatmap
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Write the interactive HTML spectrogram heatmap.
 *
 * The HTML is self-contained: STFT data is embedded as a 2-D JS array
 * and rendered with Plotly.js Heatmap (Viridis colorscale, −80..0 dB).
 * Open the file in any modern browser -- no local server needed.
 * -----------------------------------------------------------------------*/
static int write_html(const char *path,
                      const double *spec_db,
                      unsigned num_frames, unsigned num_bins,
                      const double *times_s, const double *freqs_hz)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    fprintf(fp,
        "<!DOCTYPE html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <title>Spectrogram – miniDSP</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n"
        "  <style>\n"
        "    * { box-sizing: border-box; margin: 0; padding: 0; }\n"
        "    body {\n"
        "      font-family: system-ui, -apple-system, 'Segoe UI', sans-serif;\n"
        "      background: #fafafa; color: #222; padding: 1.5rem;\n"
        "    }\n"
        "    h1 { font-size: 1.4rem; margin-bottom: 0.3rem; }\n"
        "    .subtitle {\n"
        "      color: #666; font-size: 0.9rem; margin-bottom: 1.2rem;\n"
        "    }\n"
        "    .info {\n"
        "      background: #f0f4ff; border-left: 4px solid #2563eb;\n"
        "      padding: 0.8rem 1rem; font-size: 0.85rem; line-height: 1.6;\n"
        "      border-radius: 0 6px 6px 0; max-width: 800px; margin-top: 1rem;\n"
        "    }\n"
        "    .info code { background: #e2e8f0; padding: 1px 5px; border-radius: 3px; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n"
        "  <h1>Spectrogram</h1>\n"
        "  <p class=\"subtitle\">\n"
        "    Linear chirp 200 &rarr; 4000 Hz &nbsp;|&nbsp;\n"
        "    2 s at 16 kHz &nbsp;|&nbsp;\n"
        "    FFT N=512 (32 ms) &nbsp;|&nbsp;\n"
        "    hop=128 (75%% overlap)\n"
        "  </p>\n"
        "  <div id=\"spectrogram\"></div>\n"
        "  <div class=\"info\">\n"
        "    <strong>How to read this plot:</strong><br>\n"
        "    Each column is one STFT frame (one windowed FFT). The colour shows\n"
        "    magnitude in dB: <code>20 log10(|X|/N)</code>, floored at &minus;80 dB.\n"
        "    The diagonal stripe is the chirp sweeping from 200 Hz at t=0 to\n"
        "    4000 Hz at t=2 s. Hover over the heatmap to read exact values.\n"
        "  </div>\n\n");

    /* Embed time axis */
    fprintf(fp, "  <script>\n");
    fprintf(fp, "    const times = [");
    for (unsigned f = 0; f < num_frames; f++) {
        fprintf(fp, "%.4f", times_s[f]);
        if (f + 1 < num_frames) fprintf(fp, ",");
    }
    fprintf(fp, "];\n");

    /* Embed frequency axis */
    fprintf(fp, "    const freqs = [");
    for (unsigned k = 0; k < num_bins; k++) {
        fprintf(fp, "%.2f", freqs_hz[k]);
        if (k + 1 < num_bins) fprintf(fp, ",");
    }
    fprintf(fp, "];\n");

    /* Embed spectrogram as 2-D array: z[k][f] for Plotly heatmap
     * (Plotly heatmap: z is indexed as z[row][col] where row = y, col = x) */
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
        "      hovertemplate: 't: %%{x:.3f} s<br>f: %%{y:.0f} Hz<br>%%{z:.1f} dB<extra></extra>'\n"
        "    }], {\n"
        "      xaxis: { title: 'Time (s)' },\n"
        "      yaxis: { title: 'Frequency (Hz)' },\n"
        "      margin: { t: 30, r: 80, b: 50, l: 70 },\n"
        "      height: 500\n"
        "    }, { responsive: true });\n"
        "  </script>\n"
        "</body>\n"
        "</html>\n");

    fclose(fp);
    return 0;
}

int main(void)
{
    /* ------------------------------------------------------------------
     * Signal parameters
     * ----------------------------------------------------------------*/
    const double   sample_rate  = 16000.0;   /* Hz                        */
    const double   duration_s   = 2.0;       /* seconds                   */
    const unsigned signal_len   = (unsigned)(sample_rate * duration_s); /* 32000 */
    const double   chirp_f0     = 200.0;     /* start frequency (Hz)      */
    const double   chirp_f1     = 4000.0;    /* end frequency (Hz)        */

    /* STFT parameters */
    const unsigned N   = 512;  /* window size: 32 ms at 16 kHz     */
    const unsigned hop = 128;  /* hop size: 8 ms (75% overlap)     */

    const unsigned num_bins   = N / 2 + 1;
    const unsigned num_frames = MD_stft_num_frames(signal_len, N, hop);

    /* ------------------------------------------------------------------
     * Allocate buffers
     * ----------------------------------------------------------------*/
    double *signal   = malloc(signal_len * sizeof(double));
    double *mag_out  = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *spec_db  = malloc((size_t)num_frames * num_bins * sizeof(double));
    double *times_s  = malloc(num_frames * sizeof(double));
    double *freqs_hz = malloc(num_bins   * sizeof(double));

    if (!signal || !mag_out || !spec_db || !times_s || !freqs_hz) {
        fprintf(stderr, "allocation failed\n");
        return 1;
    }

    /* ------------------------------------------------------------------
     * Generate a linear chirp: x(t) = sin(2*pi*(f0 + (f1-f0)/(2*T)*t)*t)
     *
     * A linear chirp sweeps instantaneous frequency from f0 to f1
     * linearly over duration T.  It is the canonical test signal for
     * spectrograms because the diagonal stripe it produces is immediately
     * recognisable and easy to verify visually.
     * ----------------------------------------------------------------*/
    //! [generate-chirp]
    const double chirp_rate = (chirp_f1 - chirp_f0) / duration_s;  /* Hz/s */
    for (unsigned i = 0; i < signal_len; i++) {
        double t = (double)i / sample_rate;
        signal[i] = sin(2.0 * M_PI * (chirp_f0 + 0.5 * chirp_rate * t) * t);
    }
    //! [generate-chirp]

    /* ------------------------------------------------------------------
     * Compute the STFT.
     *
     * MD_stft() applies a Hanning window to each frame internally, so
     * we do not need to pre-window the signal.  The output is a row-major
     * matrix: mag_out[f * num_bins + k] = |X_f(k)| (unnormalised FFTW).
     * ----------------------------------------------------------------*/
    //! [compute-stft]
    MD_stft(signal, signal_len, N, hop, mag_out);
    //! [compute-stft]

    /* ------------------------------------------------------------------
     * Convert magnitudes to dB.
     *
     * Normalise by N before taking the log so that a full-scale sine
     * (amplitude 1) reads near 0 dB.  Floor at 1e-6 (= -120 dB) to
     * prevent log(0).
     * ----------------------------------------------------------------*/
    //! [convert-db]
    for (unsigned i = 0; i < num_frames * num_bins; i++) {
        spec_db[i] = 20.0 * log10(fmax(mag_out[i] / (double)N, 1e-6));
    }
    //! [convert-db]

    /* Build time and frequency axes */
    for (unsigned f = 0; f < num_frames; f++) {
        times_s[f] = (double)(f * hop) / sample_rate;
    }
    for (unsigned k = 0; k < num_bins; k++) {
        freqs_hz[k] = (double)k * sample_rate / (double)N;
    }

    /* ------------------------------------------------------------------
     * Write results to CSV
     * ----------------------------------------------------------------*/
    const char *csv_file = "spectrogram.csv";
    FILE *fp = fopen(csv_file, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", csv_file);
        return 1;
    }

    fprintf(fp, "frame,time_s,bin,freq_hz,magnitude_db\n");
    for (unsigned f = 0; f < num_frames; f++) {
        for (unsigned k = 0; k < num_bins; k++) {
            fprintf(fp, "%u,%.4f,%u,%.2f,%.2f\n",
                    f, times_s[f], k, freqs_hz[k],
                    spec_db[f * num_bins + k]);
        }
    }
    fclose(fp);
    printf("Wrote %u frames x %u bins to %s\n", num_frames, num_bins, csv_file);

    /* ------------------------------------------------------------------
     * Write interactive HTML heatmap
     * ----------------------------------------------------------------*/
    const char *html_file = "spectrogram.html";
    if (write_html(html_file, spec_db, num_frames, num_bins,
                   times_s, freqs_hz) == 0) {
        printf("Wrote interactive spectrogram to %s\n", html_file);
    }

    printf("Chirp: %.0f Hz -> %.0f Hz over %.1f s\n",
           chirp_f0, chirp_f1, duration_s);
    printf("STFT: N=%u, hop=%u, %u frames, %u bins\n",
           N, hop, num_frames, num_bins);

    /* ------------------------------------------------------------------
     * Cleanup
     * ----------------------------------------------------------------*/
    free(freqs_hz);
    free(times_s);
    free(spec_db);
    free(mag_out);
    free(signal);
    MD_shutdown();

    return 0;
}
