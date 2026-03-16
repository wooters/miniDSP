## 1. Update dB floor in C source files

- [x] 1.1 Change `zmin` from `-80` to `-100` in spectext spectrogram calls in `docs/gen_signal_plots.c`
- [x] 1.2 Change `zmin` from `-80` to `-100` in HTML output in `examples/spectrogram_text.c`
- [x] 1.3 Change PNG dB-to-index mapping floor from `-80` to `-100` in `examples/spectrogram_text.c`

## 2. Regenerate docs assets

- [x] 2.1 Rebuild `gen_signal_plots` and run to regenerate `guides/plots/spectext_hello_spectrogram.html` and `guides/plots/steg_spectext_spectrogram.html`

## 3. Verify

- [x] 3.1 Confirm regenerated HTML files contain `zmin: -100`
- [x] 3.2 Build and run `spectrogram_text` to verify HTML and PNG output use -100 dB floor
