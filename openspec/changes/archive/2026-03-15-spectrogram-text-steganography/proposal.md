## Why

miniDSP's audio steganography module currently offers two methods: LSB (high-capacity, fragile) and frequency-band BFSK (robust, low-capacity). Both are machine-readable only — there's no way to visually verify that a message is present without running `decode()`. A hybrid method that combines LSB data encoding with spectrogram text art in the ultrasonic band (18-24 kHz) gives users the best of both worlds: reliable machine-readable round-trip via LSB, plus a human-readable visual watermark visible in any spectrogram viewer. The visual channel also serves as a quick tamper indicator — if the spectrogram text is intact, the LSB data likely is too.

The library already has the two building blocks needed: `MD_spectrogram_text()` for bitmap font rendering into audio, and `MD_resample()` for sample rate conversion. This change composes them into a new steganography method.

## What Changes

- Add `MD_STEG_SPECTEXT` (value `2`) to the steganography method constants
- Implement spectext encode path in `src/minidsp_steg.c`: upsample host to 48 kHz if needed, LSB-encode the message, generate spectrogram text art at 18-24 kHz with fixed 30 ms column width, mix at low amplitude (~0.02), write 48 kHz output
- Implement spectext decode path: delegate to existing `lsb_decode()` (the LSB data channel carries the payload)
- Implement spectext detection in `MD_steg_detect()`: LSB header present + ultrasonic energy distinguishes from plain LSB
- Update `MD_steg_capacity()` for the new method: capacity limited by whichever channel is smaller (visual capacity based on host duration and fixed column width, or LSB capacity based on sample count)
- Update CLI tool `tools/audio_steg/audio_steg.c` to accept `"spectext"` as a method name
- Add tests for round-trip encode/decode, capacity calculation, resampling trigger, visual truncation, and detection
- Expand the audio steganography guide with a new section covering the spectext method, including formulas, "Reading the formula in C" snippets, listening comparison, and spectrogram visualization
- Generate new multimedia assets: stego WAV file, spectrogram plot showing the embedded text, and listening comparison audio

## Capabilities

### New Capabilities
- `spectext-steg`: Hybrid LSB + spectrogram text art steganography method with automatic upsampling to 48 kHz

### Modified Capabilities
- `steg-detection`: `MD_steg_detect()` updated to distinguish `MD_STEG_SPECTEXT` from `MD_STEG_LSB` by probing for ultrasonic energy
- `steg-cli`: `audio_steg` CLI tool accepts `"spectext"` method and handles 48 kHz output

## Impact

- **Core library**: `src/minidsp_steg.c` gains spectext encode/decode/detect/capacity paths; new `#include` dependency on spectrogram text and resampler functions (already in `libminidsp.a`)
- **Public API**: One new constant `MD_STEG_SPECTEXT`; no new function signatures — existing `MD_steg_encode`, `MD_steg_decode`, `MD_steg_capacity`, `MD_steg_detect` all gain a new method branch
- **CLI tool**: `tools/audio_steg/audio_steg.c` updated for new method name and self-test
- **Tests**: New test functions in `tests/test_steg.c`
- **Documentation**: Expanded `guides/audio-steganography.md` with new section; new audio/plot assets in `guides/audio/` and `guides/plots/`; `Doxyfile` `HTML_EXTRA_FILES` updated
- **Output sample rate**: Unlike LSB/BFSK which preserve the host sample rate, spectext always outputs 48 kHz WAV (upsamples if necessary)
