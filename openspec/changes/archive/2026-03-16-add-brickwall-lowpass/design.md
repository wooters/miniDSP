## Context

The spectext steganography encoder upsamples host audio to 48 kHz using a polyphase sinc resampler (`MD_resample()` with 32 zero-crossings, Kaiser beta=10). This resampler achieves ~100 dB stopband attenuation but has a ~1.5 kHz transition bandwidth starting at the original Nyquist. For a 16 kHz host (Nyquist 8 kHz), spectral images leak through from 8-9.5 kHz and propagate faintly to higher bands, competing with the spectext tones at 18-23.5 kHz.

The library already has FFTW3 as a dependency and uses r2c/c2r FFT plans in `minidsp_spectrum.c` (r2c only) and `minidsp_gcc.c` (both r2c and c2r).

## Goals / Non-Goals

**Goals:**
- Add a general-purpose brickwall lowpass filter function with zero transition bandwidth
- Eliminate resampler spectral images in the spectext encode pipeline
- Make spectext messages clearly readable in any spectrogram viewer without special display settings

**Non-Goals:**
- Changing the resampler itself (increasing zero-crossings, adding filter options)
- Adding other filter types (highpass, bandpass, Butterworth, etc.)
- Changing spectrogram visualization parameters (window type, dB range) — the fix is in the signal, not the viewer
- Cached FFT plan infrastructure for this function (one-off plans are sufficient)

## Decisions

### 1. FFT-based brickwall over higher-order FIR

**Decision**: Implement as FFT → zero bins above cutoff → IFFT.

**Alternatives considered**:
- *More resampler zero-crossings*: Narrows the transition band but never eliminates it. Doesn't solve the fundamental problem.
- *High-order FIR lowpass*: Essentially duplicates what the resampler already does, just with more taps. More code for a worse result.
- *Biquad cascade*: Limited stopband attenuation without many stages. Not practical for brickwall behavior.

**Rationale**: FFT brickwall gives mathematically perfect suppression (zero energy above cutoff, zero transition bandwidth) with ~15 lines of code. The only downside — Gibbs ringing at the cutoff — is irrelevant here because the cutoff (e.g., 8 kHz for 16 kHz input) is far from the spectext band (18 kHz).

### 2. In-place operation

**Decision**: `MD_lowpass_brickwall(signal, len, cutoff_hz, sample_rate)` modifies the buffer directly.

**Alternatives considered**:
- *Separate input/output buffers*: More flexible but unnecessary — the use case is always "clean this buffer before the next step." Extra allocation for no benefit.

**Rationale**: Matches the pipeline use case. The caller (spectext_encode) already has a mutable buffer (`mixed[]`). In-place avoids an extra malloc/free cycle.

### 3. One-off FFT plans (not cached)

**Decision**: Create and destroy r2c/c2r plans within each call to `MD_lowpass_brickwall()`.

**Alternatives considered**:
- *Share the spectrum module's cached r2c plan and add a paired c2r*: Would complicate the plan cache, which is keyed on `_spec_N` and currently only manages r2c. The brickwall's FFT size equals the signal length, which differs from the spectrum module's typical usage (STFT window size).
- *Add a separate brickwall plan cache*: More infrastructure for a function called once per encode.

**Rationale**: `MD_lowpass_brickwall()` is called once per spectext encode. FFTW plan creation with `FFTW_ESTIMATE` is fast (no measurement). The simplicity of one-off plans outweighs the minor cost. If performance becomes critical in future, caching can be added without API change.

### 4. Lives in minidsp_spectrum.c

**Decision**: Place the implementation in `minidsp_spectrum.c` alongside other FFT-based functions.

**Rationale**: Follows the module boundary convention — FFT-dependent APIs belong in `minidsp_spectrum.c`. The function uses `fftw_plan_dft_r2c_1d`, `fftw_plan_dft_c2r_1d`, and `fftw_execute`, all of which are already used in the spectrum module's sibling `minidsp_gcc.c`.

### 5. Cutoff bin calculation: round down

**Decision**: The cutoff bin is `floor(cutoff_hz / (sample_rate / N))`. All bins from cutoff_bin+1 through N/2 are zeroed.

**Rationale**: Conservative — ensures no energy above the cutoff frequency. The DC bin (0) and bins up to and including the cutoff bin are preserved.

## Risks / Trade-offs

**[Gibbs ringing at cutoff]** → The brickwall filter produces ringing artifacts (overshoot/undershoot) near the cutoff frequency. For the spectext use case, the cutoff is at the original Nyquist (e.g., 8 kHz) and the spectext band starts at 18 kHz — 10 kHz of separation. The ringing is inaudible and invisible in the spectext region.

**[Memory allocation for FFT buffers]** → The function allocates a complex buffer of size (N/2+1) for the FFT. For typical audio lengths (48000 samples/sec × 3 sec = 144000 samples), this is ~1.1 MB. Acceptable for a one-shot offline operation.

**[FFTW thread safety]** → `fftw_plan_*` functions are not thread-safe. This matches the existing library contract — miniDSP functions are not advertised as thread-safe. No change needed.
