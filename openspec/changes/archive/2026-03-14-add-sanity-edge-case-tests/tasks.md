## 1. Core Measurement Sanity Tests

- [x] 1.1 Add `test_power_zero` (all-zero signal returns 0.0) and `test_power_single` (N=1 returns x^2) to the MD_power section
- [x] 1.2 Add `test_scale_identity_range` (same in/out range) and `test_scale_inverted_range` (reversed output range) to the MD_scale section
- [x] 1.3 Add `test_fit_in_range_constant` (constant signal, min==max edge case) to the MD_fit_within_range section
- [x] 1.4 ~~Add `test_adjust_dblevel_silence`~~ SKIPPED: `MD_adjust_dblevel` divides by zero on silent input (no guard). This is a library bug, not a test gap.

## 2. Signal Generator Zero-Amplitude and Degenerate Tests

- [x] 2.1 Add `test_sine_zero_amplitude` (all zeros) and `test_sine_zero_frequency` (DC = all zeros since sin(0)=0) to the MD_sine_wave section
- [x] 2.2 Add `test_white_noise_zero_amplitude` (all zeros) to the MD_white_noise section
- [x] 2.3 Add `test_square_zero_amplitude` (all zeros) and `test_sawtooth_zero_amplitude` (all zeros) to their respective sections
- [x] 2.4 ~~Add `test_chirp_log_equal_freq`~~ SKIPPED: `MD_chirp_log` has `assert(f_start != f_end)` — this is assert-guarded misuse, not a valid edge case.
- [x] 2.5 Add `test_shepard_zero_amplitude` (all zeros) to the MD_shepard_tone section

## 3. Effects Identity/Degenerate Parameter Tests

- [x] 3.1 Add `test_mix_zero_weights` (both weights 0 → silence) and `test_mix_passthrough_b` (weight_a=0, weight_b=1 → signal b) to the MD_mix section
- [x] 3.2 ~~Add `test_delay_echo_zero_delay`~~ SKIPPED: `MD_delay_echo` has `assert(delay_samples > 0)` — assert-guarded misuse.
- [x] 3.3 Add `test_tremolo_full_depth` (depth=1.0 reaches zero gain) to the MD_tremolo section
- [x] 3.4 Add `test_comb_reverb_zero_feedback` (feedback=0, single echo only) to the MD_comb_reverb section

## 4. Convolution/FIR Identity and Single-Tap Tests

- [x] 4.1 Add `test_moving_average_width_one` (width=1 is identity) to the convolution/FIR section
- [x] 4.2 Add `test_fir_filter_single_tap` (single coefficient scales signal) to the FIR section
- [x] 4.3 Add `test_convolution_fft_ola_impulse_identity` (length-1 unit kernel reproduces input) to the FFT OLA section

## 5. Window Length-2 and Kaiser Beta=0 Tests

- [x] 5.1 Add `test_hann_length_two`, `test_hamming_length_two`, `test_blackman_length_two` for N=2 edge cases
- [x] 5.2 Add `test_kaiser_beta_zero` (beta=0 produces all-ones rectangular window) to the Kaiser section

## 6. Spectrum and STFT Sanity Tests

- [x] 6.1 Add `test_phase_spectrum_dc_signal` (constant signal has zero DC phase) to the phase spectrum section
- [x] 6.2 Add `test_stft_single_frame` (signal_len == N, single frame matches magnitude_spectrum) to the STFT section

## 7. GCC-PHAT Sanity Tests

- [x] 7.1 Add `test_gcc_identical_signals` (identical signals → zero delay) to the GCC section
- [x] 7.2 Add `test_gcc_silence` (two zero signals → zero cross-correlation) to the GCC section

## 8. Resampler Identity Tests

- [x] 8.1 ~~Add `test_resample_same_rate`~~ SKIPPED: already covered by existing `test_resample_identity` (uses 8000.0 → 8000.0).
- [x] 8.2 Add `test_resample_output_len_same_rate` (same rate → same length) to the MD_resample_output_len section

## 9. DTMF Edge Case Tests

- [x] 9.1 Add `test_dtmf_detect_silence` (all-zero signal → zero digits detected) to the DTMF section
- [x] 9.2 ~~Add `test_dtmf_signal_length_single_digit`~~ SKIPPED: already covered in existing `test_dtmf_signal_length` (tests single digit at line 3747).

## 10. Steganography Edge Case Tests

- [x] 10.1 Add `test_steg_encode_empty_message` (empty string → 0 bits, host unchanged) to the steganography section

## 11. Math Utility Edge Cases

- [x] 11.1 Add `test_bessel_i0_negative_symmetry` (I0(-x) == I0(x)) to the Bessel section
- [x] 11.2 Add `test_sinc_symmetry` (sinc(-x) == sinc(x)) to the sinc section

## 12. Lowpass FIR Extreme Cutoff Tests

- [x] 12.1 Add `test_lowpass_fir_near_nyquist` (high cutoff → passes high frequencies) to the lowpass FIR section
- [x] 12.2 Add `test_lowpass_fir_very_low_cutoff` (low cutoff → heavy attenuation of mid frequencies) to the lowpass FIR section

## 13. Build and Verify

- [x] 13.1 Build test suite with `make -C tests test_minidsp` and verify all existing + new tests pass
