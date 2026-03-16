## MODIFIED Requirements

### Requirement: Spectext spectrograms use Hot colormap by default
The spectext spectrogram HTML output (both `gen_signal_plots.c` docs assets and `spectrogram_text.c` example program) SHALL use the Plotly `'Hot'` colorscale instead of `'Viridis'` for spectrogram text art visualizations. The dB range SHALL be `zmin: -100, zmax: 0`.

#### Scenario: HTML spectrogram uses Hot colorscale with -100 dB floor
- **WHEN** `gen_signal_plots` generates `spectext_hello_spectrogram.html` or `steg_spectext_spectrogram.html`
- **THEN** the Plotly heatmap trace SHALL specify `colorscale: 'Hot'`
- **AND** `zmin` SHALL be `-100` and `zmax` SHALL be `0`

#### Scenario: Example program HTML output uses Hot colorscale with -100 dB floor
- **WHEN** `spectrogram_text` generates `spectrogram_text.html`
- **THEN** the Plotly heatmap trace SHALL specify `colorscale: 'Hot'`
- **AND** `zmin` SHALL be `-100` and `zmax` SHALL be `0`

### Requirement: PNG output uses Hot colormap LUT
The `spectrogram_text.c` example program SHALL use a 256-entry `Hot` colormap lookup table for PNG output instead of the Viridis LUT. The linear dB-to-index mapping logic SHALL use the range [-100, 0] dB.

#### Scenario: PNG spectrogram renders with Hot colors and -100 dB floor
- **WHEN** `spectrogram_text` generates `spectrogram_text.png` (default colormap)
- **THEN** pixels SHALL be colored using the Hot colormap (black → red → yellow → white)
- **AND** the dB-to-index mapping SHALL be linear over [-100, 0] dB

#### Scenario: Grayscale option still works
- **WHEN** `spectrogram_text --colormap grayscale` is run
- **THEN** pixels SHALL be colored using linear grayscale with the [-100, 0] dB range

### Requirement: Docs HTML assets are regenerated
The tracked HTML files `guides/plots/spectext_hello_spectrogram.html` and `guides/plots/steg_spectext_spectrogram.html` SHALL be regenerated with the new -100 dB floor after the C code changes.

#### Scenario: Regenerated spectext docs assets
- **WHEN** `make docs` is run after the code changes
- **THEN** `guides/plots/spectext_hello_spectrogram.html` SHALL contain `zmin: -100`
- **AND** `guides/plots/steg_spectext_spectrogram.html` SHALL contain `zmin: -100`
