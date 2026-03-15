## Context

miniDSP has FIR filtering (`MD_fir_filter`, `MD_convolution_*`), window functions (Hann, Hamming, Blackman, Rect), and signal generators, but no sample rate conversion. A consultant PRD proposes a polyphase sinc resampler with supporting math utilities. This design refines that PRD into an implementation plan aligned with existing codebase conventions.

Existing infrastructure:
- `MD_dot` in `minidsp_core.c` — dot product, the inner loop of polyphase filtering
- `MD_Gen_*_Win` family — window generators following the `(double *out, unsigned n)` pattern
- `MD_fir_filter` / `MD_convolution_*` in `minidsp_fir.c` — FIR filtering infrastructure
- No persistent state needed — offline resampling allocates/frees per call

## Goals / Non-Goals

**Goals:**
- High-quality offline sample rate conversion (>100 dB SNR with default parameters)
- Support arbitrary rate pairs (44100↔48000, 48000↔16000, 44100↔22050, etc.)
- Composable building blocks: Bessel I₀, sinc, Kaiser window, and FIR lowpass design are independently useful public functions
- Follow existing API patterns: caller-allocated output buffers, `assert()` preconditions, `MD_` prefix

**Non-Goals:**
- Streaming/block-based API with state carry-over (future work)
- Multi-channel support (caller loops over channels)
- SIMD optimization (correctness first)
- Exact rational-ratio optimization via GCD decomposition (unnecessary complexity for offline use — a fixed-phase approach with interpolation handles all cases uniformly)

## Decisions

### 1. Unified fixed-phase interpolation instead of GCD-based rational/irrational split

The PRD proposes GCD reduction for rational ratios and fixed phases for irrational ratios. This design uses a single approach: a fixed number of sub-phases (512) with linear interpolation between adjacent phases for all rate ratios.

**Why:** Simpler implementation, same quality for offline use, follows "clarity over cleverness." The memory difference is trivial (512 × 64 taps = 32K doubles = 256 KB for 32 zero crossings). Exact rational optimization can be added later if profiling shows a need.

### 2. Polyphase filter table: internal allocation with malloc/free

`MD_resample` allocates a `num_phases × (2 * num_zero_crossings)` filter table internally, uses it, and frees it before returning. No caching across calls.

**Why:** Rate ratios change between calls, so caching provides no benefit (unlike FFT plans which are reused for same-length signals). This avoids adding teardown complexity to `MD_shutdown()`.

**Alternative considered:** Caller-allocated filter table. Rejected because it exposes internal structure and the allocation is fast relative to the filtering work.

### 3. MD_design_lowpass_fir is independent of the resampler

The resampler builds its own polyphase table by evaluating the Kaiser-windowed sinc at specific fractional offsets per phase. It does NOT call `MD_design_lowpass_fir`. Both are public but serve different use cases.

**Why:** `MD_design_lowpass_fir` produces a single contiguous FIR filter (useful for anti-aliasing, band-limiting). The resampler needs per-phase evaluation at fractional sample offsets, which is structurally different.

### 4. File placement follows module boundaries

| Function | File | Rationale |
|---|---|---|
| `MD_bessel_i0`, `MD_sinc`, `MD_Gen_Kaiser_Win` | `src/minidsp_core.c` | Stateless math utilities |
| `MD_design_lowpass_fir` | `src/minidsp_fir.c` | FIR design alongside existing FIR functions |
| `MD_resample`, `MD_resample_output_len` | `src/minidsp_resample.c` (new) | New module for resampling domain |

### 5. Anti-aliasing cutoff and gain correction

- Cutoff frequency: `min(in_rate, out_rate) / 2` — prevents aliasing when downsampling and limits bandwidth when upsampling
- For downsampling, the sinc filter coefficients are scaled by `out_rate / in_rate` to preserve signal energy (the downsampling rate factor)
- Boundary handling: clamp input indices to `[0, input_len-1]` conceptually zero-padding beyond boundaries

### 6. Fixed 512 sub-phases with linear interpolation

For each output sample:
1. Compute fractional input position: `pos = n * in_rate / out_rate`
2. Integer part: `idx = floor(pos)` — center input sample
3. Fractional part: `frac = pos - idx` — selects between sub-phases
4. Phase index: `phase = floor(frac * num_phases)`, interpolation weight: `alpha = frac * num_phases - phase`
5. Interpolate between `table[phase]` and `table[phase+1]` coefficients
6. Dot product with `2 * num_zero_crossings` input samples centered around `idx`

512 phases gives ~60 dB interpolation accuracy, well above the Kaiser window's stopband attenuation for typical β values.

## Risks / Trade-offs

**[Memory allocation per call]** → Acceptable for offline use. The table is small (~256 KB worst case) and allocated once. If streaming API is added later, the table can be cached in a state struct.

**[No SIMD]** → Performance is O(output_len × 2 × num_zero_crossings) multiplications. For typical audio (1 second at 48 kHz, 64 taps per phase), that's ~6M multiplications — fast enough without SIMD.

**[Fixed 512 phases may be excessive for simple ratios]** → Negligible overhead. Table construction is O(512 × 2 × num_zero_crossings) = O(32K) operations, dwarfed by the filtering loop.

**[Edge effects at signal boundaries]** → Zero-padding at boundaries introduces a brief fade-in/fade-out. This is standard behavior for offline resamplers. Documented in the API.
