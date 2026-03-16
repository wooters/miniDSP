## 1. Library Function: MD_lowpass_brickwall()

- [x] 1.1 Add `MD_lowpass_brickwall()` declaration and Doxygen doc-comment to `include/minidsp.h` (brief, formula, params, code example)
- [x] 1.2 Update the feature list at the top of `include/minidsp.h` to include FFT-based brickwall lowpass filtering
- [x] 1.3 Implement `MD_lowpass_brickwall()` in `src/minidsp_spectrum.c`: assert preconditions, allocate complex buffer, create one-off r2c and c2r plans (FFTW_ESTIMATE), execute r2c, zero bins above cutoff, execute c2r, normalize by N, free plans and buffer

## 2. Integrate into Spectext Encode

- [x] 2.1 In `src/minidsp_steg.c` `spectext_encode()`, increase `MD_resample()` zero-crossings from 32 to 128 to narrow the transition bandwidth
- [x] 2.2 Add a call to `MD_lowpass_brickwall(mixed, out_len, sample_rate / 2.0, SPECTEXT_TARGET_SR)` after the upsampling block and before the spectrogram text mixing — only when upsampling was performed (`sample_rate < SPECTEXT_TARGET_SR`)

## 3. Tests

- [x] 3.1 Add test `test_lowpass_brickwall_preserves_low_freq` to `tests/test_spectrum.c`: generate a 1000 Hz sine at 48 kHz, apply brickwall at 8000 Hz, verify energy is preserved within 0.1%
- [x] 3.2 Add test `test_lowpass_brickwall_removes_high_freq`: generate a 20000 Hz sine at 48 kHz, apply brickwall at 8000 Hz, verify energy is effectively zero
- [x] 3.3 Add test `test_lowpass_brickwall_mixed_signal`: generate a signal with 1000 Hz + 20000 Hz tones, apply brickwall at 8000 Hz, verify only the 1000 Hz component remains
- [x] 3.4 Register the new tests with `RUN_TEST()` in `run_spectrum_tests()` under a `printf("\n--- MD_lowpass_brickwall ---\n")` section header
- [x] 3.5 Run `make test` and verify all tests pass (existing and new)

## 4. Documentation

- [x] 4.1 Update the encode pipeline diagram in `guides/audio-steganography.md` to show `MD_lowpass_brickwall()` between upsampling and spectrogram text mixing
- [x] 4.2 Regenerate the spectext spectrogram asset by running `make docs` and verify the text is clearly readable

## 5. Verification

- [x] 5.1 Run full test suite (`make test`) to confirm no regressions
- [x] 5.2 Build docs (`make docs`) and visually inspect the regenerated spectext spectrogram
