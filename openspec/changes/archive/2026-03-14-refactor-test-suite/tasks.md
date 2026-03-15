## 1. Shared test infrastructure

- [x] 1.1 Create `tests/test_helpers.h` with `extern` counter declarations, `RUN_TEST` macro, `approx_equal` static inline, and `delay_signal` static helper
- [x] 1.2 Create a forward-declaration header or convention: each test file declares its `run_*_tests()` prototype (driver includes them)

## 2. Extract per-module test files

- [x] 2.1 Create `tests/test_core.c` — extract tests for dot, energy, power, power_db, scale, fit_within_range, adjust_dblevel, entropy, rms, zcr, autocorrelation, peak_detect, pitch, mix. Expose `run_core_tests(void)`
- [x] 2.2 Create `tests/test_effects.c` — extract tests for delay_echo, tremolo, comb_reverb. Expose `run_effects_tests(void)`
- [x] 2.3 Create `tests/test_generators.c` — extract tests for sine, white_noise, impulse, chirp_linear, chirp_log, square, sawtooth, shepard_tone. Expose `run_generators_tests(void)`
- [x] 2.4 Create `tests/test_spectrum.c` — extract tests for magnitude_spectrum, PSD, phase_spectrum, STFT, mel/MFCC (including shared `spectrum_fn_t` helpers). Expose `run_spectrum_tests(void)`
- [x] 2.5 Create `tests/test_fir.c` — extract tests for convolution, FIR filter, moving average, FFT OLA, bessel_i0, sinc, Kaiser window, lowpass FIR design. Expose `run_fir_tests(void)`
- [x] 2.6 Create `tests/test_gcc.c` — extract GCC-PHAT tests. Expose `run_gcc_tests(void)`
- [x] 2.7 Create `tests/test_biquad.c` — extract biquad filter tests. Expose `run_biquad_tests(void)`
- [x] 2.8 Create `tests/test_windows.c` — extract window generation tests (Hann, Hamming, Blackman, Rect). Expose `run_windows_tests(void)`
- [x] 2.9 Create `tests/test_dtmf.c` — extract DTMF tests. Expose `run_dtmf_tests(void)`
- [x] 2.10 Create `tests/test_spectext.c` — extract spectrogram text tests. Expose `run_spectext_tests(void)`
- [x] 2.11 Create `tests/test_steg.c` — extract steganography tests. Expose `run_steg_tests(void)`
- [x] 2.12 Create `tests/test_fileio.c` — extract file I/O writer tests. Expose `run_fileio_tests(void)`
- [x] 2.13 Create `tests/test_resample.c` — extract resampling tests. Expose `run_resample_tests(void)`

## 3. Rewrite driver and update build

- [x] 3.1 Rewrite `tests/test_minidsp.c` as slim driver: define global counters, call all `run_*_tests()`, print summary, call `MD_shutdown()`
- [x] 3.2 Update `tests/Makefile` to compile all `test_*.c` files to `.o` and link into `test_minidsp`

## 4. Verify

- [x] 4.1 Run `make -C tests clean && make -C tests test` — confirm 276/276 pass
- [x] 4.2 Verify output order: all 276 tests present, grouped by module (section reordering is expected per design)
