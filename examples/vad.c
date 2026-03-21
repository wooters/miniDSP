/**
 * @file vad.c
 * @brief Example: voice activity detection on a synthesized signal.
 *
 * This program demonstrates the VAD module by:
 *   1. Synthesizing a signal with alternating speech (sine bursts) and silence.
 *   2. Initializing VAD state with default params.
 *   3. Optionally calibrating on known silence frames.
 *   4. Processing frame-by-frame, printing per-frame results to CSV.
 *   5. Generating an interactive HTML visualization (Plotly.js).
 *
 * Build and run (from the repository root):
 *   make                          # build libminidsp.a first
 *   make -C examples plot         # build example, run it, generate HTML
 *   open examples/vad_plot.html    # interactive VAD timeline
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minidsp.h"
#include "plot_html.h"

/* Signal parameters */
#define SAMPLE_RATE  16000.0
#define DURATION     2.0
#define FRAME_SIZE   256
#define SPEECH_FREQ  1000.0
#define SPEECH_AMP   0.5

/* -----------------------------------------------------------------------
 * Synthesize a test signal: silence, speech, silence, speech, silence
 * -----------------------------------------------------------------------*/
static void synthesize_signal(double *signal, unsigned total_samples)
{
    /* Fill with silence first */
    for (unsigned i = 0; i < total_samples; i++)
        signal[i] = 0.0;

    /* Speech segments: [0.2s–0.6s] and [1.0s–1.5s] */
    unsigned seg1_start = (unsigned)(0.2 * SAMPLE_RATE);
    unsigned seg1_end   = (unsigned)(0.6 * SAMPLE_RATE);
    unsigned seg2_start = (unsigned)(1.0 * SAMPLE_RATE);
    unsigned seg2_end   = (unsigned)(1.5 * SAMPLE_RATE);

    for (unsigned i = seg1_start; i < seg1_end && i < total_samples; i++)
        signal[i] = SPEECH_AMP * sin(2.0 * M_PI * SPEECH_FREQ * i / SAMPLE_RATE);

    for (unsigned i = seg2_start; i < seg2_end && i < total_samples; i++)
        signal[i] = SPEECH_AMP * sin(2.0 * M_PI * SPEECH_FREQ * i / SAMPLE_RATE);
}

/* -----------------------------------------------------------------------
 * Write the interactive HTML visualization
 * -----------------------------------------------------------------------*/
static int write_html(const char *path,
                      const double *times, const double *waveform,
                      const int *decisions, const double *scores,
                      const double *feat_matrix,
                      unsigned num_frames, double threshold)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "cannot open %s for writing\n", path);
        return -1;
    }

    char subtitle[256];
    snprintf(subtitle, sizeof(subtitle),
        "%.0f Hz sample rate &nbsp;|&nbsp;\n"
        "    frame size %d &nbsp;|&nbsp;\n"
        "    duration %.1f s",
        SAMPLE_RATE, FRAME_SIZE, DURATION);

    plot_html_begin(fp, "Voice Activity Detection", subtitle, 0);

    /* Embed data arrays */
    fprintf(fp, "  <div id=\"waveform\" class=\"plot-container\"></div>\n");
    fprintf(fp, "  <div id=\"features\" class=\"plot-container\"></div>\n");
    fprintf(fp, "  <div id=\"score\" class=\"plot-container\"></div>\n");
    fprintf(fp, "  <div id=\"decision\" class=\"plot-container\"></div>\n\n");

    fprintf(fp, "  <script>\n");

    plot_html_js_array(fp, "times", times, num_frames, "%.4f");

    /* Waveform: use first sample of each frame */
    plot_html_js_array(fp, "waveform", waveform, num_frames, "%.6f");
    plot_html_js_array(fp, "scores", scores, num_frames, "%.6f");

    /* Decisions as doubles for plotting */
    fprintf(fp, "    const decisions = [");
    for (unsigned i = 0; i < num_frames; i++) {
        fprintf(fp, "%d", decisions[i]);
        if (i + 1 < num_frames) fprintf(fp, ",");
    }
    fprintf(fp, "];\n");

    /* Feature traces */
    const char *feat_names[] = {
        "energy", "zcr", "spectral_entropy",
        "spectral_flatness", "band_energy_ratio"
    };
    for (int f = 0; f < MD_VAD_NUM_FEATURES; f++) {
        double feat_vals[num_frames];
        for (unsigned i = 0; i < num_frames; i++)
            feat_vals[i] = feat_matrix[i * MD_VAD_NUM_FEATURES + f];
        plot_html_js_array(fp, feat_names[f], feat_vals, num_frames, "%.6f");
    }

    /* Plotly traces */
    fprintf(fp,
        "\n"
        "    // Waveform plot\n"
        "    Plotly.newPlot('waveform', [{\n"
        "      x: times, y: waveform, type: 'scatter', mode: 'lines',\n"
        "      line: {color: '#2563eb', width: 1}, name: 'Waveform'\n"
        "    }], {\n"
        "      title: 'Waveform', height: 200,\n"
        "      margin: {l:50, r:20, t:40, b:30},\n"
        "      xaxis: {title: ''}, yaxis: {title: 'Amplitude'}\n"
        "    }, {responsive: true});\n\n"

        "    // Feature traces\n"
        "    const feat_colors = ['#2563eb','#dc2626','#16a34a','#9333ea','#ea580c'];\n"
        "    const feat_labels = ['Energy','ZCR','Spectral Entropy',\n"
        "                         'Spectral Flatness','Band Energy Ratio'];\n"
        "    const feat_data = [energy, zcr, spectral_entropy,\n"
        "                       spectral_flatness, band_energy_ratio];\n"
        "    const feat_traces = feat_data.map((d,i) => ({\n"
        "      x: times, y: d, type: 'scatter', mode: 'lines',\n"
        "      line: {color: feat_colors[i], width: 1.5},\n"
        "      name: feat_labels[i]\n"
        "    }));\n"
        "    Plotly.newPlot('features', feat_traces, {\n"
        "      title: 'Normalized Features (0-1)', height: 250,\n"
        "      margin: {l:50, r:20, t:40, b:30},\n"
        "      xaxis: {title: ''}, yaxis: {title: 'Value', range: [-0.05, 1.05]}\n"
        "    }, {responsive: true});\n\n"

        "    // Score with threshold\n"
        "    Plotly.newPlot('score', [\n"
        "      {x: times, y: scores, type: 'scatter', mode: 'lines',\n"
        "       line: {color: '#2563eb', width: 2}, name: 'Score'},\n"
        "      {x: [times[0], times[times.length-1]],\n"
        "       y: [%.4f, %.4f],\n"
        "       type: 'scatter', mode: 'lines',\n"
        "       line: {color: '#dc2626', width: 1.5, dash: 'dash'},\n"
        "       name: 'Threshold'}\n"
        "    ], {\n"
        "      title: 'Combined Score vs Threshold', height: 200,\n"
        "      margin: {l:50, r:20, t:40, b:30},\n"
        "      xaxis: {title: ''}, yaxis: {title: 'Score', range: [-0.05, 1.05]}\n"
        "    }, {responsive: true});\n\n"

        "    // Decision timeline\n"
        "    Plotly.newPlot('decision', [{\n"
        "      x: times, y: decisions, type: 'scatter', mode: 'lines',\n"
        "      fill: 'tozeroy', line: {color: '#16a34a', width: 1.5},\n"
        "      fillcolor: 'rgba(22,163,74,0.3)', name: 'Decision'\n"
        "    }], {\n"
        "      title: 'Speech Decision (1=speech, 0=silence)', height: 150,\n"
        "      margin: {l:50, r:20, t:40, b:40},\n"
        "      xaxis: {title: 'Time (s)'},\n"
        "      yaxis: {title: '', range: [-0.1, 1.2], tickvals: [0,1], ticktext: ['Silence','Speech']}\n"
        "    }, {responsive: true});\n",
        threshold, threshold);

    fprintf(fp, "  </script>\n");
    plot_html_end(fp);

    fclose(fp);
    printf("  HTML: %s\n", path);
    return 0;
}

/* -----------------------------------------------------------------------
 * Main
 * -----------------------------------------------------------------------*/
int main(void)
{
    unsigned total_samples = (unsigned)(DURATION * SAMPLE_RATE);
    double *signal = malloc(total_samples * sizeof(double));
    synthesize_signal(signal, total_samples);

    unsigned num_frames = total_samples / FRAME_SIZE;

    /* Allocate output arrays */
    double *times     = malloc(num_frames * sizeof(double));
    double *waveform  = malloc(num_frames * sizeof(double));
    int    *decisions = malloc(num_frames * sizeof(int));
    double *scores    = malloc(num_frames * sizeof(double));
    double *feat_matrix = malloc(num_frames * MD_VAD_NUM_FEATURES * sizeof(double));

    /* --- Initialize VAD --- */
    //! [vad-init]
    MD_vad_params params;
    MD_vad_default_params(&params);
    params.threshold = 0.3;

    MD_vad_state state;
    MD_vad_init(&state, &params);
    //! [vad-init]

    /* --- Calibrate on first few frames (known silence) --- */
    //! [vad-calibrate]
    unsigned cal_frames = 10;
    for (unsigned i = 0; i < cal_frames; i++) {
        double *frame = signal + i * FRAME_SIZE;
        MD_vad_calibrate(&state, frame, FRAME_SIZE, SAMPLE_RATE);
    }
    //! [vad-calibrate]

    /* --- Process frame-by-frame --- */
    //! [vad-process]
    FILE *csv = fopen("vad.csv", "w");
    fprintf(csv, "frame,time,decision,score,energy,zcr,"
                 "spectral_entropy,spectral_flatness,band_energy_ratio\n");

    for (unsigned i = 0; i < num_frames; i++) {
        double *frame = signal + i * FRAME_SIZE;
        double score;
        double features[MD_VAD_NUM_FEATURES];

        int decision = MD_vad_process_frame(&state, frame, FRAME_SIZE,
                                            SAMPLE_RATE, &score, features);

        double t = (double)(i * FRAME_SIZE) / SAMPLE_RATE;
        times[i] = t;
        waveform[i] = frame[0];
        decisions[i] = decision;
        scores[i] = score;
        for (int f = 0; f < MD_VAD_NUM_FEATURES; f++)
            feat_matrix[i * MD_VAD_NUM_FEATURES + f] = features[f];

        fprintf(csv, "%u,%.4f,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                i, t, decision, score,
                features[0], features[1], features[2],
                features[3], features[4]);
    }
    fclose(csv);
    printf("  CSV: vad.csv (%u frames)\n", num_frames);
    //! [vad-process]

    /* --- Demonstrate custom weights --- */
    //! [vad-custom-weights]
    MD_vad_params energy_only_params;
    MD_vad_default_params(&energy_only_params);
    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++)
        energy_only_params.weights[i] = 0.0;
    energy_only_params.weights[MD_VAD_FEAT_ENERGY] = 1.0;

    MD_vad_state energy_only_state;
    MD_vad_init(&energy_only_state, &energy_only_params);
    //! [vad-custom-weights]

    /* Write HTML visualization */
    write_html("vad_plot.html", times, waveform, decisions, scores,
               feat_matrix, num_frames, params.threshold);

    /* Clean up */
    free(signal);
    free(times);
    free(waveform);
    free(decisions);
    free(scores);
    free(feat_matrix);
    MD_shutdown();

    return 0;
}
