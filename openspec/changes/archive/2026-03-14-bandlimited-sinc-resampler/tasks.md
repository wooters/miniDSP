## 1. Math Utilities (minidsp_core.c)

- [x] 1.1 Implement `MD_bessel_i0` — power series I₀(x) with 1e-15 relative convergence
- [x] 1.2 Implement `MD_sinc` — normalized sinc with 1e-12 near-zero threshold
- [x] 1.3 Add tests for `MD_bessel_i0`: I₀(0)=1, known values at x=1 and x=5, monotonic increase
- [x] 1.4 Add tests for `MD_sinc`: sinc(0)=1, integer zeros, sinc(0.5)=2/π, near-zero threshold

## 2. Kaiser Window (minidsp_core.c)

- [x] 2.1 Implement `MD_Gen_Kaiser_Win` — Kaiser window using I₀, assert n>0, handle n==1
- [x] 2.2 Add tests: single-sample degenerate case, symmetry, tapered ends, peak at center, β comparison

## 3. Lowpass FIR Design (minidsp_fir.c)

- [x] 3.1 Implement `MD_design_lowpass_fir` — windowed-sinc with Kaiser, normalize to unity DC gain
- [x] 3.2 Add tests: coefficient sum=1.0, linear phase symmetry, passband/stopband behavior

## 4. Resampler Core (src/minidsp_resample.c — new file)

- [x] 4.1 Implement `MD_resample_output_len` — ceil(input_len * out_rate / in_rate)
- [x] 4.2 Implement `MD_resample` — build 512-phase polyphase table, interpolate between phases, clamp boundaries
- [x] 4.3 Add tests for `MD_resample_output_len`: upsample, downsample, non-integer ratio cases
- [x] 4.4 Add tests for `MD_resample`: identity (1:1), DC preservation, sine frequency preservation
- [x] 4.5 Add tests for `MD_resample`: output length matches prediction for common rate pairs
- [x] 4.6 Add tests for `MD_resample`: energy preservation (white noise RMS within 0.5 dB)
- [x] 4.7 Add tests for `MD_resample`: anti-aliasing suppression (20 kHz sine at 48k → 16k)

## 5. Header Declarations (include/minidsp.h)

- [x] 5.1 Add Doxygen declarations for `MD_bessel_i0`, `MD_sinc`, `MD_Gen_Kaiser_Win` in the core/window sections (with `@param`, `@return`, `@code` examples, LaTeX formulas)
- [x] 5.2 Add Doxygen declaration for `MD_design_lowpass_fir` in the FIR section
- [x] 5.3 Add Doxygen declarations for `MD_resample` and `MD_resample_output_len` in a new resampling section
- [x] 5.4 Update `@brief` feature list (~lines 5-17) — add Kaiser window, sinc, Bessel I₀, lowpass FIR design, and sample rate conversion

## 6. Build System

- [x] 6.1 Add `src/minidsp_resample.c` to `MD_SRCS` in root Makefile
- [x] 6.2 Verify full build (`make clean && make`) and test suite passes (`make -C tests`)

## 7. Documentation — Guides

- [x] 7.1 Add Kaiser window section to existing `guides/window-functions.md` (formula, "Reading the formula in C:", API, `\snippet`, visuals)
- [x] 7.2 Add lowpass FIR design section to existing `guides/fir-convolution.md` (formula, "Reading the formula in C:", API, `\snippet`, visuals)
- [x] 7.3 Create new `guides/resampling.md` covering Bessel I₀, sinc, and the polyphase resampler (formulas, "Reading the formula in C:", API, `\snippet`, visuals)
- [x] 7.4 Add `\subpage resampling` entry to `guides/tutorials.md`

## 8. Documentation — Examples

- [x] 8.1 Add Kaiser window snippet to existing `examples/window_functions.c` (with `//! [snippet-name]` markers)
- [x] 8.2 Add lowpass FIR design snippet to existing `examples/fir_convolution.c` (with `//! [snippet-name]` markers)
- [x] 8.3 Create new `examples/resampler.c` — demonstrate resampling between common rate pairs (with snippet markers for guide embedding)
- [x] 8.4 Add `resampler` to `examples/Makefile` `EXAMPLES` variable and `plot:` target
- [x] 8.5 Add `.gitignore` and `.dockerignore` entries for `examples/resampler`, `examples/resampler.csv`, `examples/resampler.html`

## 9. Documentation — Doxygen Config

- [x] 9.1 Add any new HTML assets (iframe plots) to `Doxyfile` `HTML_EXTRA_FILES` — N/A (no generated HTML plots for this feature)

## 10. Finalize

- [x] 10.1 Update `README.md` — add resampling to feature list / API summary
- [x] 10.2 Bump VERSION patch number (0.1.0 → 0.2.0)
