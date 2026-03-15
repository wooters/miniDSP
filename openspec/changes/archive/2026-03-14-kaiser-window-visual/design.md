## Context

The docs guide `guides/window-functions.md` covers five window functions: Hanning, Hamming, Blackman, Rectangular, and Kaiser. The first four each have a pair of interactive Plotly.js iframes (time-domain taps + frequency-domain magnitude spectrum). Kaiser has the math, API docs, and a snippet example but no visual plots.

The plot generation infrastructure already exists — `examples/gen_signal_plots.c` has `write_window_time_html()` and `write_window_spectrum_html()` helper functions that produce iframe-ready HTML. The other four windows use these helpers with shared constants (`WINDOW_N = 256`, `WINDOW_FFT_VIS = 4096`).

## Goals / Non-Goals

**Goals:**
- Add Kaiser window time-domain and spectrum plots to the docs, visually consistent with the other four windows
- Follow the existing generation, embedding, and build pipeline exactly

**Non-Goals:**
- Adding multiple beta values or interactive beta sliders (single representative beta is sufficient)
- Adding Kaiser window tests (separate concern)
- Modifying the Kaiser library implementation

## Decisions

**Beta value: 10.0** — This matches the example snippet already in the guide (`MD_Gen_Kaiser_Win(kaiser, N, 10.0)`) and the `window_functions.c` example program. Beta=10 produces a visually distinctive shape (narrow main lobe, strong sidelobe suppression) that contrasts well with the other windows in the guide.

**Reuse existing helper functions** — `write_window_time_html()` and `write_window_spectrum_html()` already handle the Plotly layout, axis labels, and iframe-friendly CSS. No new helper code needed.

**Same parameters as other windows** — Use `WINDOW_N` (256) and `WINDOW_FFT_VIS` (4096) for consistency in visual comparison.

## Risks / Trade-offs

No meaningful risks — this is additive documentation with zero library code changes.
