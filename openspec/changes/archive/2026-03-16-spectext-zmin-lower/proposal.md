## Why

The spectext spectrogram dB floor is currently set to -80 dB. Lowering it to -100 dB will reveal more detail in the quieter portions of the spectrogram, making the text art patterns more visible against the background.

## What Changes

- Change `zmin` from `-80` to `-100` in all spectext spectrogram HTML outputs (both `gen_signal_plots.c` docs assets and `spectrogram_text.c` example program).
- Change the PNG dB-to-index mapping range from `[-80, 0]` to `[-100, 0]` in `spectrogram_text.c`.
- Regenerate tracked HTML docs assets with the new dB floor.

## Capabilities

### New Capabilities

(none)

### Modified Capabilities

- `spectext-colormap`: Change `zmin` from `-80` to `-100` for spectext spectrogram visualizations (both HTML and PNG).

## Impact

- `docs/gen_signal_plots.c` — HTML spectrogram `zmin` values
- `examples/spectrogram_text.c` — HTML `zmin` and PNG dB mapping floor
- `guides/plots/spectext_hello_spectrogram.html` — regenerated asset
- `guides/plots/steg_spectext_spectrogram.html` — regenerated asset
