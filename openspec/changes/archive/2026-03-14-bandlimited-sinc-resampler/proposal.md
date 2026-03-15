## Why

Sample rate conversion is a fundamental DSP operation missing from miniDSP. A polyphase sinc resampler is the gold-standard approach used by libsamplerate and professional audio tools — it achieves >100 dB SNR with reasonable filter lengths, handles arbitrary rate ratios, and requires no external dependencies. Adding it rounds out the library's signal processing toolkit alongside existing FIR filtering, FFT analysis, and signal generation.

## What Changes

- Add `MD_bessel_i0` — zeroth-order modified Bessel function I₀(x), a general-purpose math utility
- Add `MD_sinc` — normalized sinc function sin(πx)/(πx)
- Add `MD_Gen_Kaiser_Win` — Kaiser window generator (joins existing Hann/Hamming/Blackman/Rect family)
- Add `MD_design_lowpass_fir` — Kaiser-windowed sinc FIR lowpass filter design
- Add `MD_resample` — polyphase sinc resampler for offline sample rate conversion
- Add `MD_resample_output_len` — helper to compute required output buffer size
- New source file `src/minidsp_resample.c` for the resampler core
- Build system updated to include the new source file

## Capabilities

### New Capabilities
- `math-utilities`: Bessel I₀ and normalized sinc — standalone math functions useful beyond resampling
- `kaiser-window`: Kaiser window generator with configurable β parameter for sidelobe control
- `lowpass-fir-design`: Windowed-sinc FIR lowpass filter design with Kaiser window and configurable cutoff
- `polyphase-resampler`: Offline polyphase sinc resampler supporting arbitrary rate ratios

### Modified Capabilities

(none)

## Impact

- **New source file**: `src/minidsp_resample.c` added to `MD_SRCS` in root Makefile
- **Public API**: 6 new functions declared in `include/minidsp.h`
- **Existing modules modified**: `src/minidsp_core.c` (3 new functions), `src/minidsp_fir.c` (1 new function)
- **Tests**: New test functions in `tests/test_minidsp.c` for all 6 functions
- **Dependencies**: None new — uses only standard C math library
- **Memory**: `MD_resample` allocates/frees a polyphase filter table internally per call (no persistent state, no changes to `MD_shutdown`)
