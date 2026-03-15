## ADDED Requirements

### Requirement: Core measurement functions have zero/single-element sanity tests
Every core measurement function (`MD_power`, `MD_power_db`, `MD_scale`, `MD_scale_vec`, `MD_fit_within_range`, `MD_adjust_dblevel`) SHALL have tests for trivial inputs where the expected output is self-evident.

#### Scenario: MD_power of zero signal returns zero
- **WHEN** `MD_power` is called with an all-zero signal
- **THEN** it SHALL return 0.0

#### Scenario: MD_power of single element returns x^2
- **WHEN** `MD_power` is called with N=1 and signal `{5.0}`
- **THEN** it SHALL return 25.0

#### Scenario: MD_scale identity range returns input unchanged
- **WHEN** `MD_scale` is called with in_min==out_min and in_max==out_max
- **THEN** the output SHALL equal the input

#### Scenario: MD_scale inverted output range reverses direction
- **WHEN** `MD_scale` is called mapping [0,10] to [100,0]
- **THEN** input 0 SHALL map to 100 and input 10 SHALL map to 0

#### Scenario: MD_fit_within_range constant signal maps to range midpoint
- **WHEN** `MD_fit_within_range` is called with a constant signal (all same values)
- **THEN** all output values SHALL be the midpoint of the target range (since min==max, range collapses)

#### Scenario: MD_adjust_dblevel with silent input produces silent output
- **WHEN** `MD_adjust_dblevel` is called with an all-zero signal
- **THEN** the output SHALL remain all zeros (no amplification of silence)

### Requirement: Signal generators handle zero-amplitude and degenerate parameters
Signal generator functions SHALL produce correct output for zero-amplitude and other degenerate-but-valid parameter combinations.

#### Scenario: MD_sine_wave with zero amplitude produces all zeros
- **WHEN** `MD_sine_wave` is called with amplitude=0.0
- **THEN** all output samples SHALL be 0.0

#### Scenario: MD_sine_wave with zero frequency produces DC signal
- **WHEN** `MD_sine_wave` is called with freq=0.0 and amplitude=1.0
- **THEN** all output samples SHALL be 0.0 (sin(0)=0 at all samples)

#### Scenario: MD_white_noise with zero amplitude produces all zeros
- **WHEN** `MD_white_noise` is called with amplitude=0.0
- **THEN** all output samples SHALL be 0.0

#### Scenario: MD_square_wave with zero amplitude produces all zeros
- **WHEN** `MD_square_wave` is called with amplitude=0.0
- **THEN** all output samples SHALL be 0.0

#### Scenario: MD_sawtooth_wave with zero amplitude produces all zeros
- **WHEN** `MD_sawtooth_wave` is called with amplitude=0.0
- **THEN** all output samples SHALL be 0.0

#### Scenario: MD_chirp_log with equal start and end frequency behaves like constant-frequency
- **WHEN** `MD_chirp_log` is called with f_start == f_end
- **THEN** the output SHALL have energy concentrated at that single frequency

#### Scenario: MD_shepard_tone with zero amplitude produces all zeros
- **WHEN** `MD_shepard_tone` is called with amplitude=0.0
- **THEN** all output samples SHALL be 0.0

### Requirement: Effects functions handle identity/degenerate parameters
Effect functions (`MD_delay_echo`, `MD_tremolo`, `MD_comb_reverb`, `MD_mix`) SHALL produce correct output when parameters reduce the effect to an identity or null operation.

#### Scenario: MD_mix with both weights zero produces silence
- **WHEN** `MD_mix` is called with weight_a=0.0 and weight_b=0.0
- **THEN** all output samples SHALL be 0.0

#### Scenario: MD_mix with weight_b=1.0 passes through second signal
- **WHEN** `MD_mix` is called with weight_a=0.0 and weight_b=1.0
- **THEN** output SHALL equal signal b exactly

#### Scenario: MD_delay_echo with delay=0 produces wet+dry passthrough
- **WHEN** `MD_delay_echo` is called with delay_samples=0
- **THEN** output SHALL equal dry_gain*input + wet_gain*input (immediate echo)

#### Scenario: MD_tremolo at full depth reaches zero gain
- **WHEN** `MD_tremolo` is called with depth=1.0
- **THEN** the gain envelope SHALL reach 0.0 at its minimum

#### Scenario: MD_comb_reverb with feedback=0 produces single echo
- **WHEN** `MD_comb_reverb` is called with feedback=0.0 on an impulse
- **THEN** the output SHALL contain the impulse and exactly one echo (no further repeats)

### Requirement: Convolution and FIR functions handle identity kernels and single-tap cases
Convolution functions SHALL produce correct output for identity kernels (single unit impulse) and width-1/single-tap degenerate cases.

#### Scenario: MD_moving_average with width=1 is identity
- **WHEN** `MD_moving_average` is called with width=1
- **THEN** output SHALL equal the input signal exactly

#### Scenario: MD_fir_filter with single tap scales signal
- **WHEN** `MD_fir_filter` is called with a single coefficient `{0.5}`
- **THEN** output SHALL be input scaled by 0.5

#### Scenario: MD_convolution_fft_ola with length-1 impulse kernel reproduces input
- **WHEN** `MD_convolution_fft_ola` is called with kernel `{1.0}` (length 1)
- **THEN** output SHALL equal the input signal

### Requirement: Window functions handle length=2 edge case
All window generation functions SHALL produce valid output for the minimum useful length N=2.

#### Scenario: Hann window of length 2 has both endpoints zero
- **WHEN** `MD_Gen_Hann_Win` is called with n=2
- **THEN** both output values SHALL be 0.0

#### Scenario: Hamming window of length 2 has both endpoints equal to alpha0-alpha1
- **WHEN** `MD_Gen_Hamming_Win` is called with n=2
- **THEN** both endpoints SHALL equal the Hamming minimum value (~0.08)

#### Scenario: Blackman window of length 2 has both endpoints near zero
- **WHEN** `MD_Gen_Blackman_Win` is called with n=2
- **THEN** both endpoints SHALL be near 0.0

#### Scenario: Kaiser window with beta=0 equals rectangular window
- **WHEN** `MD_Gen_Kaiser_Win` is called with beta=0.0
- **THEN** all output values SHALL be 1.0 (since I0(0)=1)

### Requirement: Spectrum functions handle constant/DC signals as sanity checks
FFT-based spectrum functions SHALL produce correct, self-evident output for constant (DC-only) signals where the spectrum is trivially known.

#### Scenario: MD_phase_spectrum of constant signal has zero phase at DC
- **WHEN** `MD_phase_spectrum` is called with a constant signal
- **THEN** the DC bin (index 0) phase SHALL be 0.0 (or very near it)

#### Scenario: MD_stft single-frame case matches magnitude spectrum
- **WHEN** `MD_stft` is called with signal_len == N (one frame, no overlap)
- **THEN** the output SHALL match `MD_magnitude_spectrum` applied to the same windowed frame

### Requirement: GCC-PHAT handles identical and silent signal pairs
`MD_get_delay` and `MD_gcc` SHALL produce correct results for trivial signal pairs.

#### Scenario: MD_get_delay with identical signals returns zero delay
- **WHEN** `MD_get_delay` is called with two copies of the same signal
- **THEN** the estimated delay SHALL be 0

#### Scenario: MD_gcc with silent signals produces zero cross-correlation
- **WHEN** `MD_gcc` is called with two all-zero signals
- **THEN** all cross-correlation output values SHALL be 0.0

### Requirement: Resampler handles identity rate and single-sample input
`MD_resample` SHALL produce correct output for trivial resampling cases.

#### Scenario: MD_resample at same rate is near-identity
- **WHEN** `MD_resample` is called with input_rate == output_rate
- **THEN** output length SHALL equal input length and values SHALL closely match the input

#### Scenario: MD_resample_output_len with same rate returns same length
- **WHEN** `MD_resample_output_len` is called with input_rate == output_rate
- **THEN** the returned length SHALL equal input_len

### Requirement: DTMF handles edge cases
DTMF functions SHALL handle degenerate-but-valid inputs.

#### Scenario: MD_dtmf_detect on silent signal returns zero digits
- **WHEN** `MD_dtmf_detect` is called with an all-zero signal
- **THEN** the number of detected digits SHALL be 0

#### Scenario: MD_dtmf_signal_length for single digit returns correct length
- **WHEN** `MD_dtmf_signal_length` is called with num_digits=1
- **THEN** the returned length SHALL equal one tone duration plus one pause duration in samples

### Requirement: Steganography handles empty message
Steganography encode/decode functions SHALL handle zero-length messages gracefully.

#### Scenario: MD_steg_encode with empty message returns zero bits written
- **WHEN** `MD_steg_encode` is called with an empty string ("")
- **THEN** the function SHALL return 0 and output SHALL equal host signal

#### Scenario: MD_steg_capacity returns capacity for valid signal
- **WHEN** `MD_steg_capacity` is called with a signal of known length
- **THEN** the returned capacity SHALL be positive and proportional to signal length

### Requirement: Math utility edge cases
`MD_bessel_i0` and `MD_sinc` SHALL handle extreme inputs correctly.

#### Scenario: MD_bessel_i0 with negative input returns same as positive
- **WHEN** `MD_bessel_i0` is called with a negative value
- **THEN** the result SHALL equal `MD_bessel_i0(|x|)` (I0 is an even function)

#### Scenario: MD_sinc is symmetric
- **WHEN** `MD_sinc` is called with x and -x
- **THEN** both results SHALL be equal (sinc is an even function)

### Requirement: Lowpass FIR design handles extreme cutoffs
`MD_design_lowpass_fir` SHALL handle cutoff frequencies at extremes.

#### Scenario: MD_design_lowpass_fir near-Nyquist cutoff approaches all-pass
- **WHEN** `MD_design_lowpass_fir` is called with cutoff near Nyquist
- **THEN** the filter SHALL pass a high-frequency tone with minimal attenuation

#### Scenario: MD_design_lowpass_fir very low cutoff attenuates most frequencies
- **WHEN** `MD_design_lowpass_fir` is called with a very low cutoff (e.g., 100 Hz at 44100 Hz)
- **THEN** a 1000 Hz tone filtered through it SHALL be heavily attenuated
