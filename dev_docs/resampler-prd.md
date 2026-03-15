# PRD: Bandlimited Sinc Resampler for miniDSP

## Overview

Add high-quality sample rate conversion to miniDSP using a polyphase Kaiser-windowed sinc interpolation algorithm. This requires several new utility functions (Bessel I₀, Kaiser window, normalized sinc, FIR lowpass design) plus the resampler core.

## Motivation

Sample rate conversion is a fundamental DSP operation missing from the library. A polyphase sinc resampler is the gold-standard approach (used by libsamplerate and professional audio tools), achieves >100 dB SNR with reasonable filter lengths, handles arbitrary input/output rate ratios, and requires no external dependencies beyond what miniDSP already uses.

## Codebase Context

### Existing infrastructure to build on

- **Window functions** in `src/minidsp_core.c`: Hann, Hamming, Blackman, Rectangular (`MD_Gen_*_Win`)
- **FIR filtering** in `src/minidsp_fir.c`: `MD_fir_filter`, `MD_convolution_time`, `MD_convolution_fft_ola`
- **`MD_dot`** in `src/minidsp_core.c`: dot product — this is the inner loop of a polyphase FIR
- **FFT infrastructure** via FFTW with cached plans (`MD_shutdown` cleanup pattern)
- **Signal generators** for test signal creation (`MD_sine_wave`, `MD_chirp_linear`, etc.)

### Conventions (from CLAUDE.md)

- C17 standard (`-std=c17 -Wall -Wextra -pedantic`)
- Use `assert()` for structural misuse; sentinel returns for valid-but-unresolved runtime outcomes
- Module boundaries: time-domain/stateless analysis → `src/minidsp_core.c`; FFT-dependent → `src/minidsp_spectrum.c`
- All public API uses the `MD_` prefix
- All functions declared in `include/minidsp.h` with full Doxygen documentation including `@param`, `@return`, `@code` examples, and LaTeX formulas where appropriate
- Memory: caller-allocated output buffers (no hidden mallocs in the public API, except for internal caching patterns like the FFT plan cache)
- Header uses `#ifndef` guards, defines `M_PI` if missing

### File organization

- Source files: `src/minidsp_*.c` (one per module domain)
- Headers: `include/minidsp.h` (single public header)
- Tests: `tests/test_minidsp.c` (single test file, linked against libminidsp.a + FFTW3 + libsndfile + math)
- Build: `Makefile` at root; `MD_SRCS` list must include any new `.c` file

## New Functions

### 1. `MD_bessel_i0` — Zeroth-order modified Bessel function of the first kind

**File:** `src/minidsp_core.c` (stateless math utility)

```c
double MD_bessel_i0(double x);
```

- Compute I₀(x) using the standard power series: `I₀(x) = Σ [(x/2)^k / k!]²` for k = 0, 1, 2, ...
- Converge until the term is below `1e-15 * sum` (relative tolerance)
- This is a pure function with no allocations

**Why it's public:** I₀ is useful beyond Kaiser windows (e.g., Kaiser-Bessel derived windows, statistical distributions). Making it a named utility follows the library's composability philosophy.

### 2. `MD_Gen_Kaiser_Win` — Kaiser window generator

**File:** `src/minidsp_core.c` (alongside the other `MD_Gen_*_Win` functions)

```c
void MD_Gen_Kaiser_Win(double *out, unsigned n, double beta);
```

- Formula: `w[i] = I₀(β * sqrt(1 - ((2i/(n-1)) - 1)²)) / I₀(β)`
- For `n == 1`, output `1.0` (degenerate case)
- `beta` controls the sidelobe/mainlobe tradeoff:
  - β ≈ 5.0 → ~45 dB stopband attenuation (fast/draft quality)
  - β ≈ 10.0 → ~100 dB stopband attenuation (high quality)
  - β ≈ 14.0 → ~120 dB stopband attenuation (mastering quality)
- Assert: `n > 0`

### 3. `MD_sinc` — Normalized sinc function

**File:** `src/minidsp_core.c`

```c
double MD_sinc(double x);
```

- Returns `sin(πx) / (πx)` for `x != 0`, returns `1.0` for `x == 0`
- Use a small epsilon threshold (e.g., `fabs(x) < 1e-12`) for the zero check to handle floating-point near-zero values

### 4. `MD_design_lowpass_fir` — Windowed-sinc FIR lowpass filter design

**File:** `src/minidsp_fir.c` (alongside existing FIR functions)

```c
void MD_design_lowpass_fir(double *coeffs, unsigned num_taps,
                           double cutoff_freq, double sample_rate,
                           double kaiser_beta);
```

- Generate a Kaiser-windowed sinc lowpass filter
- `cutoff_freq` is the -6 dB point in Hz
- `num_taps` is the filter length (odd is conventional but not required)
- Normalized cutoff: `fc = cutoff_freq / sample_rate`
- For each tap i: `coeffs[i] = 2 * fc * sinc(2 * fc * (i - (num_taps-1)/2.0)) * kaiser[i]`
- Normalize coefficients so they sum to 1.0 (unity DC gain)
- Assert: `num_taps > 0`, `sample_rate > 0`, `cutoff_freq > 0`, `cutoff_freq < sample_rate / 2`

**Why this is a separate public function:** FIR lowpass design is useful independently of resampling (e.g., anti-aliasing before decimation, band-limiting signals for analysis). It makes the library more composable.

### 5. `MD_resample` — Polyphase sinc resampler

**File:** `src/minidsp_resample.c` (new source file)

```c
unsigned MD_resample(const double *input, unsigned input_len,
                     double *output, unsigned max_output_len,
                     double in_rate, double out_rate,
                     unsigned num_zero_crossings, double kaiser_beta);
```

**Parameters:**

- `input` / `input_len`: source signal
- `output` / `max_output_len`: destination buffer (caller-allocated). Use `MD_resample_output_len` to size it.
- `in_rate` / `out_rate`: input and output sample rates in Hz
- `num_zero_crossings`: number of zero-crossings of the sinc on each side (controls quality). Reasonable range: 8 (fast) to 64 (high quality). Default recommendation: 32.
- `kaiser_beta`: Kaiser window β parameter. Recommendation: 10.0 for high quality.
- **Returns:** number of samples written to `output`

**Algorithm:**

1. Compute the effective rate ratio `L/M` by reducing `out_rate/in_rate` via GCD (use an internal static helper, not public API — GCD is too trivial to be a library function).
2. Compute the anti-aliasing cutoff: `fc = min(in_rate, out_rate) / 2`.
3. Build the polyphase filter table internally:
   - Total filter length: `2 * num_zero_crossings * max(L, M)` taps
   - Number of sub-phases: `L` (or a fixed resolution like 512 for irrational ratios)
   - Each sub-phase is a Kaiser-windowed sinc evaluated at the appropriate fractional offset
4. For each output sample, determine the corresponding fractional input position, select the appropriate sub-phase (or interpolate between two adjacent sub-phases for irrational ratios), and compute the dot product with the input samples.
5. Handle boundaries by zero-padding the input conceptually (no actual allocation needed — just clamp indices).

**Design decisions:**

- The function allocates the polyphase table internally and frees it before returning. This follows the library's pattern of no persistent hidden state (unlike the FFT cache, there's no performance benefit to caching resampler tables since rate ratios change between calls).
- For rational ratios (e.g., 44100→48000), use exact integer polyphase decomposition. For irrational or high-ratio cases, use a fixed table size (e.g., 512 phases) with linear interpolation between adjacent phases.
- Assert: `input_len > 0`, `in_rate > 0`, `out_rate > 0`, `num_zero_crossings > 0`, `max_output_len >= MD_resample_output_len(input_len, in_rate, out_rate)`

### 6. `MD_resample_output_len` — Compute output buffer size

```c
unsigned MD_resample_output_len(unsigned input_len,
                                double in_rate, double out_rate);
```

- Returns `ceil(input_len * out_rate / in_rate)`
- Assert: `input_len > 0`, `in_rate > 0`, `out_rate > 0`

## File Changes Summary

| File | Action |
|---|---|
| `include/minidsp.h` | Add declarations for all 6 new functions, with full Doxygen docs, in appropriate sections |
| `src/minidsp_core.c` | Add `MD_bessel_i0`, `MD_sinc`, `MD_Gen_Kaiser_Win` |
| `src/minidsp_fir.c` | Add `MD_design_lowpass_fir` |
| `src/minidsp_resample.c` | New file: `MD_resample`, `MD_resample_output_len`, internal GCD helper |
| `Makefile` | Add `src/minidsp_resample.c` to `MD_SRCS` |
| `tests/test_minidsp.c` | Add tests (see below) |

## Test Plan

All tests go in `tests/test_minidsp.c`, following the existing test patterns in that file.

### Unit tests for building blocks

1. **`MD_bessel_i0`**: verify `I₀(0) == 1.0`; spot-check against known values (e.g., `I₀(1) ≈ 1.2660658...`, `I₀(5) ≈ 27.239871...`); verify monotonic increase for x > 0
2. **`MD_sinc`**: verify `sinc(0) == 1.0`; verify `sinc(n) == 0.0` for integer n ≠ 0; verify `sinc(0.5) ≈ 2/π`
3. **`MD_Gen_Kaiser_Win`**: verify symmetry (`w[i] == w[n-1-i]`); verify `w[0] > 0` and `w[0] < w[n/2]` (tapered ends); verify peak at center; verify `n == 1` → `w[0] == 1.0`
4. **`MD_design_lowpass_fir`**: verify coefficients sum to 1.0 (DC gain); verify symmetry (linear phase); verify passband/stopband behavior by convolving with a known signal

### Integration tests for the resampler

5. **Identity**: resample at 1:1 ratio → output matches input within floating-point tolerance
6. **DC preservation**: resample a DC signal (all ones) → output is all ones
7. **Sine preservation**: generate a 440 Hz sine at 44100 Hz, resample to 48000 Hz, measure F0 of the output → should be 440 Hz within a tight tolerance (< 1 Hz)
8. **Length correctness**: verify `MD_resample_output_len` matches actual samples written for several rate pairs (44100→48000, 48000→44100, 16000→8000, 8000→16000, 44100→22050)
9. **Energy preservation**: resample white noise, compare RMS of input vs output → should be approximately equal (within 0.5 dB)
10. **Anti-aliasing**: generate a signal with content above the target Nyquist (e.g., 20 kHz sine at 48000 Hz, downsample to 16000 Hz), verify the alias is suppressed (magnitude spectrum of output should show no energy at the aliased frequency)
11. **Known rate pairs**: test the common conversions (44100↔48000, 44100↔22050, 48000↔16000) and verify output length and spectral correctness

## Out of Scope

- **Streaming/block-based API** — This PRD covers offline (full-signal) resampling only. A streaming API with state carry-over can be added later if needed.
- **Multi-channel support** — The function operates on mono signals. Multi-channel resampling is just a loop over channels at the call site.
- **SIMD optimization** — Get correctness first. The polyphase dot product is naturally vectorizable and can be optimized later with `#ifdef` blocks for SSE/NEON.
- **Inverse STFT** — Would be useful for spectral verification but is a separate feature.

## Implementation Order

1. `MD_bessel_i0` + tests
2. `MD_sinc` + tests
3. `MD_Gen_Kaiser_Win` + tests
4. `MD_design_lowpass_fir` + tests
5. `MD_resample_output_len` + tests
6. `MD_resample` + tests
7. Update Makefile
8. Bump VERSION patch number