/**
 * @file plot_html.h
 * @brief Shared HTML boilerplate for miniDSP example visualisations.
 *
 * Provides static inline helpers that emit the common DOCTYPE, CSS,
 * CDN script tags, and closing tags used by all example programs.
 * Each example only needs to write its own plot-specific content.
 */
#ifndef PLOT_HTML_H
#define PLOT_HTML_H

#include <stdio.h>

/**
 * Write the HTML preamble: DOCTYPE through opening \<body\>, CSS, CDN scripts.
 *
 * @param fp        Open file pointer.
 * @param title     Page title (shown in \<title\> and \<h1\>).
 * @param subtitle  Descriptive line below the title (raw HTML).
 * @param use_d3    If non-zero, include the D3-array CDN script.
 */
static inline void plot_html_begin(FILE *fp, const char *title,
                                   const char *subtitle, int use_d3)
{
    fprintf(fp,
        "<!DOCTYPE html>\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "  <meta charset=\"utf-8\">\n"
        "  <title>%s</title>\n"
        "  <script src=\"https://cdn.plot.ly/plotly-2.35.2.min.js\"></script>\n",
        title);
    if (use_d3) {
        fprintf(fp,
            "  <script src=\"https://d3js.org/d3-array.v3.min.js\"></script>\n");
    }
    fprintf(fp,
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
        "    .plot-container { margin-bottom: 1.5rem; }\n"
        "    .info {\n"
        "      background: #f0f4ff; border-left: 4px solid #2563eb;\n"
        "      padding: 0.8rem 1rem; font-size: 0.85rem; line-height: 1.6;\n"
        "      border-radius: 0 6px 6px 0; max-width: 700px;\n"
        "    }\n"
        "    .info code { background: #e2e8f0; padding: 1px 5px; border-radius: 3px; }\n"
        "  </style>\n"
        "</head>\n"
        "<body>\n");
    fprintf(fp, "  <h1>%s</h1>\n", title);
    fprintf(fp, "  <p class=\"subtitle\">\n    %s\n  </p>\n", subtitle);
}

/** Write the closing HTML tags. */
static inline void plot_html_end(FILE *fp)
{
    fprintf(fp,
        "</body>\n"
        "</html>\n");
}

/**
 * Emit a JavaScript array literal: const \<name\> = [v0,v1,...];
 *
 * @param fp    Open file pointer.
 * @param name  JS variable name.
 * @param data  Array of doubles to emit.
 * @param len   Number of elements.
 * @param fmt   printf format for each value (e.g. "%.4f" or "%.8f").
 */
static inline void plot_html_js_array(FILE *fp, const char *name,
                                      const double *data, unsigned len,
                                      const char *fmt)
{
    fprintf(fp, "    const %s = [", name);
    for (unsigned i = 0; i < len; i++) {
        fprintf(fp, fmt, data[i]);
        if (i + 1 < len) fprintf(fp, ",");
    }
    fprintf(fp, "];\n");
}

#endif /* PLOT_HTML_H */
