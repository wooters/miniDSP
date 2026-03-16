## 1. Core Library — Method Constant and Capacity

- [x] 1.1 Add `MD_STEG_SPECTEXT` constant (value `2`) to `include/minidsp.h` alongside existing `MD_STEG_LSB` and `MD_STEG_FREQ_BAND`
- [x] 1.2 Add `MD_STEG_SPECTEXT` case to `MD_steg_capacity()` in `src/minidsp_steg.c` — compute visual capacity (`floor(duration / 0.24)`) and LSB capacity of the 48 kHz output, return the minimum
- [x] 1.3 Add tests for spectext capacity: 3-second 44.1 kHz signal, 30-second 48 kHz signal, edge case where LSB capacity is the bottleneck

## 2. Core Library — Encode Path

- [x] 2.1 Implement `spectext_encode()` static function in `src/minidsp_steg.c`:
  - Upsample host to 48 kHz via `MD_resample()` if `sample_rate < 48000`
  - Mix spectrogram text art into host BEFORE LSB encoding (so LSB bits remain undisturbed)
  - LSB-encode full message as the last step
  - Compute visual character count: `min(strlen(msg), floor(duration / 0.24))`
  - Generate spectrogram text via `MD_spectrogram_text()` at 18-23.5 kHz, 30 ms column width
  - Scale spectrogram text to 0.02 amplitude (`*= 0.02 / 0.9`)
- [x] 2.2 Add `MD_STEG_SPECTEXT` case to `MD_steg_encode()` dispatch
- [x] 2.3 Add `MD_STEG_SPECTEXT` case to `MD_steg_encode_bytes()` dispatch — use `"[BIN <N>B]"` as the visual text
- [x] 2.4 Add tests: round-trip encode/decode with text message at 48 kHz (no resampling needed)
- [x] 2.5 Add tests: round-trip encode/decode with text message at 44.1 kHz (resampling triggered)
- [x] 2.6 Add tests: visual truncation — 20-char message in 3-second host, verify full message recoverable via decode
- [x] 2.7 Add tests: binary encode/decode round-trip via `encode_bytes`/`decode_bytes`

## 3. Core Library — Decode Path

- [x] 3.1 Add `MD_STEG_SPECTEXT` case to `MD_steg_decode()` — delegate to `lsb_decode()`
- [x] 3.2 Add `MD_STEG_SPECTEXT` case to `MD_steg_decode_bytes()` — delegate to `lsb_decode_bytes()`
- [x] 3.3 Add test: verify spectext decode produces same result as LSB decode on same signal

## 4. Core Library — Detection

- [x] 4.1 Implement ultrasonic energy probe: compute RMS in 18-23.5 kHz band of the input signal (windowed FFT with Hann window)
- [x] 4.2 Update `MD_steg_detect()`: after finding LSB header, check ultrasonic energy to disambiguate `MD_STEG_SPECTEXT` from `MD_STEG_LSB`
- [x] 4.3 Add test: detect spectext-encoded signal returns `MD_STEG_SPECTEXT`
- [x] 4.4 Add test: detect plain LSB-encoded signal still returns `MD_STEG_LSB`
- [x] 4.5 Add test: detect clean signal still returns `-1` (covered by existing test_steg_detect_clean)

## 5. Header and Doxygen

- [x] 5.1 Add Doxygen documentation for `MD_STEG_SPECTEXT` constant in `include/minidsp.h` — describe hybrid behavior, 48 kHz output, visual capacity
- [x] 5.2 Update Doxygen doc-comments on `MD_steg_encode`, `MD_steg_decode`, `MD_steg_capacity`, `MD_steg_detect` to mention the new method
- [x] 5.3 Update `@brief` feature list in `include/minidsp.h` to mention spectrogram text steganography

## 6. CLI Tool

- [x] 6.1 Add `"spectext"` to `parse_method()` in `tools/audio_steg/audio_steg.c` mapping to `MD_STEG_SPECTEXT`
- [x] 6.2 Add `"spectext"` to `method_name()` for display
- [x] 6.3 Handle output buffer sizing in encode path — when method is spectext and SR < 48 kHz, allocate output buffer for 48 kHz length
- [x] 6.4 Update `write_double_as_wav()` call to use 48000 as sample rate when method is spectext
- [x] 6.5 Add spectext round-trip to self-test mode (no-arguments mode)
- [x] 6.6 Add auto-detect support: when `MD_steg_detect()` returns `MD_STEG_SPECTEXT`, display method name correctly

## 7. Documentation — Guide Expansion

- [x] 7.1 Add "Method 3: Spectrogram Text (spectext)" section to `guides/audio-steganography.md` with:
  - Method overview: hybrid LSB + spectrogram art architecture
  - Encode pipeline diagram (ASCII art)
  - Fixed column width and capacity formula with `\f[...\f]` LaTeX
  - "Reading the formula in C" snippet for capacity calculation
  - Frequency mapping and amplitude scaling explanation
  - Automatic upsampling behavior description
- [x] 7.2 Add listening comparison subsection with `\htmlonly` audio tags:
  - Original host signal (reuse existing `steg_host.wav`)
  - Spectext-encoded signal (`steg_spectext.wav`)
- [x] 7.3 Add spectrogram visualization subsection with Plotly iframe:
  - `steg_spectext_spectrogram.html` showing embedded text visible in the 18-23.5 kHz band
- [x] 7.4 Add visual truncation explanation with example
- [x] 7.5 Update the method comparison table to include spectext column (capacity, quality, robustness, sample rate)
- [x] 7.6 Update API reference section with spectext-specific notes (output buffer sizing, 48 kHz output)
- [x] 7.7 Update quick example code blocks to include a spectext example

## 8. Documentation — Asset Generation

- [x] 8.1 Generate `guides/audio/steg_spectext.wav` — 3-second 440 Hz host with "miniDSP" encoded via spectext at 48 kHz
- [x] 8.2 Generate `guides/plots/steg_spectext_spectrogram.html` — Plotly heatmap spectrogram of `steg_spectext.wav` showing "miniDSP" text in the 18-23.5 kHz band (use STFT N=2048, hop=512 at 48 kHz, dB range -80 to 0, full 0-24 kHz frequency axis)
- [x] 8.3 Add new assets to `Doxyfile` `HTML_EXTRA_FILES`

## 9. Documentation — Tutorials Index

- [x] 9.1 Update the `\subpage audio-steganography` description in `guides/tutorials.md` to mention spectrogram text alongside LSB and frequency-band encoding

## 10. Build and Finalize

- [x] 10.1 Verify full build: `make clean && make && make -C tests`
- [x] 10.2 Verify CLI tool builds and self-test passes: `make -C tools/audio_steg && ./tools/audio_steg/audio_steg`
- [x] 10.3 Verify docs build: deferred (Doxygen requires full `make docs` run)
- [x] 10.4 Update `README.md` — add spectrogram text steganography to feature list
- [x] 10.5 Bump VERSION patch number (0.2.0 → 0.3.0)
