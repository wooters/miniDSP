## Why

When the spectext steganography method upsamples a host audio file (e.g., TIMIT at 16 kHz to 48 kHz), the polyphase sinc resampler's finite transition bandwidth (~1.5 kHz) leaves spectral images partially visible in the spectrogram. These appear as faint "mirrors" of the original content at multiples of the original Nyquist (8 kHz, 16 kHz), making the spectext tones at 18-23.5 kHz difficult to read against the image leakage. The spectext message should be clearly readable in any spectrogram viewer without special display settings.

## What Changes

- **New public API function `MD_lowpass_brickwall()`** — FFT-based brickwall lowpass filter that zeroes all frequency bins above a cutoff frequency. General-purpose utility, not spectext-specific.
- **Spectext encode pipeline updated** — increase resampler zero-crossings from 32 to 128 (quarters the transition bandwidth), then apply `MD_lowpass_brickwall()` at the original Nyquist before mixing spectext tones. The tighter resampler reduces image energy; the brickwall eliminates it completely.
- **Documentation updated** — encode pipeline diagram in the audio steganography guide reflects the new lowpass step. New function gets Doxygen doc-comment and tests.

## Capabilities

### New Capabilities
- `brickwall-lowpass`: FFT-based brickwall lowpass filter function (`MD_lowpass_brickwall`) — in-place filtering, zeroes frequency bins above cutoff, normalizes after IFFT.

### Modified Capabilities
- `spectext-steg`: Spectext encode pipeline adds a post-resampling brickwall lowpass step to eliminate spectral images before mixing spectext tones.

## Impact

- **Code**: `src/minidsp_spectrum.c` (new function + FFTW c2r plan), `src/minidsp_steg.c` (call lowpass after resampling), `include/minidsp.h` (declaration + feature list)
- **Tests**: New tests in `tests/test_minidsp_spectrum.c`
- **Docs**: `guides/audio-steganography.md` (pipeline diagram), `examples/gen_signal_plots.c` (regenerate spectext spectrogram asset)
- **APIs**: New public function. No breaking changes — existing encode/decode behavior is preserved; output quality improves.
- **Dependencies**: No new dependencies (uses existing FFTW3)
