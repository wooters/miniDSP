## ADDED Requirements

### Requirement: Kaiser window time-domain plot generation
`gen_signal_plots` SHALL generate `guides/plots/kaiser_window_time.html` containing an interactive Plotly.js chart of Kaiser window taps (beta=10.0, N=256), using the same `write_window_time_html()` helper and parameters as the other window plots.

#### Scenario: Time-domain plot file is produced
- **WHEN** `examples/gen_signal_plots` runs as part of `make docs`
- **THEN** `guides/plots/kaiser_window_time.html` is created with a Plotly chart showing 256 Kaiser window tap values

### Requirement: Kaiser window spectrum plot generation
`gen_signal_plots` SHALL generate `guides/plots/kaiser_window_spectrum.html` containing an interactive Plotly.js magnitude spectrum chart (zero-padded to 4096, dB scale), using the same `write_window_spectrum_html()` helper and parameters as the other window plots.

#### Scenario: Spectrum plot file is produced
- **WHEN** `examples/gen_signal_plots` runs as part of `make docs`
- **THEN** `guides/plots/kaiser_window_spectrum.html` is created with a Plotly chart showing the Kaiser window magnitude response in dB

### Requirement: Kaiser plots listed in Doxyfile HTML_EXTRA_FILES
Both `guides/plots/kaiser_window_time.html` and `guides/plots/kaiser_window_spectrum.html` SHALL be listed in the Doxyfile `HTML_EXTRA_FILES` so Doxygen copies them into `docs/html/`.

#### Scenario: Doxygen copies Kaiser plot files
- **WHEN** `doxygen Doxyfile` runs
- **THEN** both Kaiser plot HTML files are present in the `docs/html/` output directory

### Requirement: Kaiser guide section includes visual iframes
The Kaiser window section in `guides/window-functions.md` SHALL include `\htmlonly` iframe embeds for both the time-domain and spectrum plots, matching the pattern used by the Hanning, Hamming, Blackman, and Rectangular sections.

#### Scenario: Kaiser section displays plots in rendered docs
- **WHEN** a user views the Window Functions guide in the Doxygen HTML output
- **THEN** the Kaiser section shows two interactive iframe plots (time-domain taps and frequency spectrum)
