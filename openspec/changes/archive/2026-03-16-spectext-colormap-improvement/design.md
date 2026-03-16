## Context

All spectrogram visualizations in the miniDSP docs use Plotly's built-in `'Viridis'` colorscale with a linear dB mapping from -80 to 0. For most spectrograms (DTMF, chirp, etc.) this works well. However, for spectrogram text art — where the goal is to read embedded text — Viridis's perceptual uniformity works against us: the tones at higher frequencies (which may be slightly quieter due to the 1/N normalization spreading energy) appear as muted greens/teals rather than bright yellows, making the text hard to distinguish from the mid-range noise floor.

Two code paths produce spectext spectrograms:
1. **Plotly HTML** — `gen_signal_plots.c` (generates `spectext_hello_spectrogram.html` and `steg_spectext_spectrogram.html`) and `spectrogram_text.c` (example program)
2. **PNG** — `spectrogram_text.c` uses an embedded 256-entry Viridis LUT with linear dB-to-index mapping

## Goals / Non-Goals

**Goals:**
- Make the embedded text clearly readable in spectext spectrogram visualizations
- Change only colormap-related properties (colorscale name, dB-to-color mapping)
- Apply changes only to spectext spectrograms (not DTMF, chirp, or other docs spectrograms)

**Non-Goals:**
- Changing STFT parameters (FFT size, hop, window)
- Changing audio synthesis (frequencies, amplitudes, duration)
- Changing the Viridis colormap for non-spectext spectrograms
- Adding user-selectable colormaps at runtime (beyond the existing `--colormap` flag)

## Decisions

### 1. Use `Hot` colorscale for spectext Plotly spectrograms

**Choice:** Replace `colorscale: 'Viridis'` with `colorscale: 'Hot'` in spectext spectrogram HTML output.

**Rationale:** `Hot` maps silence (low dB) to black and tones (high dB) through red/orange to bright white/yellow. This produces high contrast between the text (bright) and the background (black), making the text immediately readable. It also evokes the classic "spectral waterfall" aesthetic that readers associate with spectrograms.

**Alternatives considered:**
- **`Inferno`** — Similar dark-to-bright mapping but tones appear as dull yellows rather than bright whites. Less contrast than `Hot`.
- **`Plasma`** — Purple-to-yellow. Better than Viridis for text but still lower contrast than `Hot` at the bright end.
- **Raise `zmin`** (e.g., -40 instead of -80) — Would compress the color range and make tones brighter, but would also wash out the background and lose detail in the noise floor. Rejected because it changes the dB mapping semantics.
- **Custom Plotly colorscale with gamma** — More control but adds complexity. `Hot` achieves the goal with zero custom code.

### 2. Replace the C-side Viridis LUT with a `Hot` LUT for PNG output

**Choice:** Replace the 256-entry Viridis RGB lookup table in `spectrogram_text.c` with a `Hot` LUT (black → red → yellow → white). Keep the linear dB-to-index mapping unchanged.

**Rationale:** PNG output should match the HTML visualization. The `Hot` colormap has a simple mathematical definition (piecewise linear in R, G, B channels), making it easy to generate programmatically or embed as a precomputed LUT.

### 3. Scope changes to spectext spectrograms only

**Choice:** Only modify the colorscale string in code paths that generate `spectext_hello_spectrogram.html`, `steg_spectext_spectrogram.html`, and the `spectrogram_text` example output. Leave all other spectrograms as Viridis.

**Rationale:** The readability problem is specific to text art where high contrast matters. Other spectrograms (DTMF, chirp, etc.) benefit from Viridis's perceptual uniformity for representing continuous spectral content.

## Risks / Trade-offs

- **Visual consistency across docs** — Spectext spectrograms will look different from other spectrograms (Viridis). This is acceptable because the purpose is fundamentally different (reading text vs. analyzing spectra).
- **PNG reconstruction inversion** — The `spectrogram_reconstruction_plan.md` describes inverting the Viridis LUT from PNG back to magnitudes. If the LUT changes, reconstruction code must match. → Mitigation: Update the dev doc to reference the `Hot` LUT. The reconstruction feature is not yet implemented, so no existing code breaks.
- **`--colormap viridis` flag** — The `spectrogram_text.c` example currently accepts `--colormap viridis|grayscale`. With `Hot` as the new default, the flag semantics change. → Mitigation: Rename the flag values to `hot|grayscale` and keep `viridis` as a deprecated alias, or simply replace the default.
