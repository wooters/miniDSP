## Why

The Kaiser window section in `guides/window-functions.md` has the math formula, API docs, and a code example — but unlike every other window (Hanning, Hamming, Blackman, Rectangular), it's missing the time-domain and spectrum iframe plots. This makes the Kaiser section feel incomplete and breaks the visual consistency of the guide.

## What Changes

- Add Kaiser window generation and HTML plot output to `examples/gen_signal_plots.c` (matching the existing `write_window_time_html` / `write_window_spectrum_html` pattern)
- Add the two new plot files to `Doxyfile` `HTML_EXTRA_FILES`
- Add `\htmlonly` iframe embeds for the Kaiser time-domain and spectrum plots in `guides/window-functions.md`

## Capabilities

### New Capabilities
- `kaiser-window-plots`: Generate and embed Kaiser window time-domain and spectrum visualizations in the docs, matching the pattern used by all other window functions

### Modified Capabilities

## Impact

- `examples/gen_signal_plots.c` — add ~5 lines to generate Kaiser window array and call the two existing plot functions
- `Doxyfile` — add two entries to `HTML_EXTRA_FILES`
- `guides/window-functions.md` — add iframe HTML blocks to the Kaiser section
- No API changes, no new dependencies, no breaking changes
