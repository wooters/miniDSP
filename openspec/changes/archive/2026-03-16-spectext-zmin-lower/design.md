## Context

The spectext spectrogram visualizations currently use a dB range of [-80, 0]. This applies to both HTML (Plotly `zmin`/`zmax`) and PNG (linear dB-to-index mapping) outputs. The value appears in `docs/gen_signal_plots.c` (docs asset generation) and `examples/spectrogram_text.c` (example program).

## Goals / Non-Goals

**Goals:**
- Lower the spectrogram dB floor from -80 to -100 across all spectext outputs
- Regenerate tracked HTML docs assets

**Non-Goals:**
- Changing non-spectext spectrograms (DTMF, chirp, resampler, brickwall, etc.)
- Changing `zmax` (stays at 0)
- Changing the Hot colormap or any other visual property

## Decisions

**Literal constant replacement**: Replace `-80` with `-100` in the specific spectext locations. No new constants or configuration needed — these are hardcoded values and should stay that way.

## Risks / Trade-offs

- [Wider dB range may make the text slightly less contrasty] → Acceptable tradeoff for seeing more spectral detail. The Hot colormap's dark-to-bright gradient handles wide ranges well.
