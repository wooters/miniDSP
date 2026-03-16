## 1. Replace Viridis LUT with Hot LUT in spectrogram_text.c

- [x] 1.1 Replace the 256-entry `viridis_lut` array in `examples/spectrogram_text.c` with a 256-entry `hot_lut` (black → red → yellow → white). The Hot colormap is piecewise linear: R ramps 0→255 over indices 0–95, G ramps 0→255 over indices 96–191, B ramps 0→255 over indices 192–255.
- [x] 1.2 Update `write_html()` in `examples/spectrogram_text.c`: change `colorscale: 'Viridis'` to `colorscale: 'Hot'`
- [x] 1.3 Rename the `use_viridis` variable/flag to `use_hot` (or a neutral name like `use_color`) and update the `--colormap` argument parsing accordingly (`hot|grayscale` instead of `viridis|grayscale`)
- [x] 1.4 Update the usage comment at the top of `spectrogram_text.c` to reflect the new `--colormap hot|grayscale` option

## 2. Update gen_signal_plots.c for spectext HTML assets

- [x] 2.1 In `examples/gen_signal_plots.c`, change `colorscale: 'Viridis'` to `colorscale: 'Hot'` in the code block that generates `spectext_hello_spectrogram.html` (around line 1396)
- [x] 2.2 In `examples/gen_signal_plots.c`, change `colorscale: 'Viridis'` to `colorscale: 'Hot'` in the code block that generates `steg_spectext_spectrogram.html` (around line 1782)
- [x] 2.3 Verify that all other spectrogram code blocks in `gen_signal_plots.c` (DTMF, chirp, resampler, brickwall, etc.) still use `colorscale: 'Viridis'`

## 3. Regenerate docs HTML assets

- [x] 3.1 Rebuild gen_signal_plots and regenerate the docs HTML assets (`make docs` or equivalent)
- [x] 3.2 Verify `guides/plots/spectext_hello_spectrogram.html` contains `colorscale: 'Hot'`
- [x] 3.3 Verify `guides/plots/steg_spectext_spectrogram.html` contains `colorscale: 'Hot'`

## 4. Visual verification

- [x] 4.1 Open the regenerated spectext spectrogram HTML files in a browser and confirm the text is clearly readable with bright colors against a dark background
- [x] 4.2 Run `cd examples && ./spectrogram_text "HELLO"` and verify the output HTML and PNG use the Hot colormap
