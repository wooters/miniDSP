### Requirement: Spectext spectrograms use Hot colormap by default
The spectext spectrogram HTML output (both `gen_signal_plots.c` docs assets and `spectrogram_text.c` example program) SHALL use the Plotly `'Hot'` colorscale instead of `'Viridis'` for spectrogram text art visualizations. The dB range (`zmin: -80, zmax: 0`) SHALL remain unchanged.

#### Scenario: HTML spectrogram uses Hot colorscale
- **WHEN** `gen_signal_plots` generates `spectext_hello_spectrogram.html` or `steg_spectext_spectrogram.html`
- **THEN** the Plotly heatmap trace SHALL specify `colorscale: 'Hot'`
- **AND** `zmin` and `zmax` SHALL remain `-80` and `0` respectively

#### Scenario: Example program HTML output uses Hot colorscale
- **WHEN** `spectrogram_text` generates `spectrogram_text.html`
- **THEN** the Plotly heatmap trace SHALL specify `colorscale: 'Hot'`

### Requirement: PNG output uses Hot colormap LUT
The `spectrogram_text.c` example program SHALL use a 256-entry `Hot` colormap lookup table for PNG output instead of the Viridis LUT. The linear dB-to-index mapping logic SHALL remain unchanged.

#### Scenario: PNG spectrogram renders with Hot colors
- **WHEN** `spectrogram_text` generates `spectrogram_text.png` (default colormap)
- **THEN** pixels SHALL be colored using the Hot colormap (black → red → yellow → white)
- **AND** the dB-to-index mapping SHALL remain linear over [-80, 0] dB

#### Scenario: Grayscale option still works
- **WHEN** `spectrogram_text --colormap grayscale` is run
- **THEN** pixels SHALL be colored using linear grayscale as before (no change)

### Requirement: Non-spectext spectrograms are unchanged
All spectrogram visualizations that are NOT spectrogram text art (DTMF, chirp, resampler, brickwall, etc.) SHALL continue using `colorscale: 'Viridis'`.

#### Scenario: DTMF spectrogram retains Viridis
- **WHEN** `gen_signal_plots` generates DTMF-related spectrogram HTML
- **THEN** the Plotly heatmap trace SHALL specify `colorscale: 'Viridis'`

### Requirement: Docs HTML assets are regenerated
The tracked HTML files `guides/plots/spectext_hello_spectrogram.html` and `guides/plots/steg_spectext_spectrogram.html` SHALL be regenerated with the new Hot colorscale after the C code changes.

#### Scenario: Regenerated spectext docs assets
- **WHEN** `make docs` is run after the code changes
- **THEN** `guides/plots/spectext_hello_spectrogram.html` SHALL contain `colorscale: 'Hot'`
- **AND** `guides/plots/steg_spectext_spectrogram.html` SHALL contain `colorscale: 'Hot'`
