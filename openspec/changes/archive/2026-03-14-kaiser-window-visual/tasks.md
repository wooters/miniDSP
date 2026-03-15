## 1. Plot Generation

- [x] 1.1 In `examples/gen_signal_plots.c`, add a Kaiser window array (`double kaiser[WINDOW_N]`), call `MD_Gen_Kaiser_Win(kaiser, WINDOW_N, 10.0)`, and pass it to `write_window_time_html()` and `write_window_spectrum_html()` to produce `guides/plots/kaiser_window_time.html` and `guides/plots/kaiser_window_spectrum.html`

## 2. Build Pipeline

- [x] 2.1 Add `guides/plots/kaiser_window_time.html` and `guides/plots/kaiser_window_spectrum.html` to the `HTML_EXTRA_FILES` list in `Doxyfile`

## 3. Guide Embedding

- [x] 3.1 Add `\htmlonly` iframe embeds for the Kaiser time-domain and spectrum plots in the Kaiser section of `guides/window-functions.md`, matching the pattern used by the other window sections
