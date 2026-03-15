## ADDED Requirements

### Requirement: Output length computation
The library SHALL provide `MD_resample_output_len(unsigned input_len, double in_rate, double out_rate)` returning `ceil(input_len * out_rate / in_rate)`.

The function SHALL assert: input_len > 0, in_rate > 0, out_rate > 0.

#### Scenario: Upsample length calculation
- **WHEN** `MD_resample_output_len(44100, 44100.0, 48000.0)` is called
- **THEN** the return value SHALL be `48000`

#### Scenario: Downsample length calculation
- **WHEN** `MD_resample_output_len(48000, 48000.0, 44100.0)` is called
- **THEN** the return value SHALL be `44100`

#### Scenario: Non-integer ratio
- **WHEN** `MD_resample_output_len(1000, 44100.0, 48000.0)` is called
- **THEN** the return value SHALL be `ceil(1000 * 48000.0 / 44100.0)` = `1089`

### Requirement: Polyphase sinc resampler
The library SHALL provide `MD_resample(input, input_len, output, max_output_len, in_rate, out_rate, num_zero_crossings, kaiser_beta)` performing offline sample rate conversion using polyphase sinc interpolation.

The function SHALL:
- Use a fixed table of 512 sub-phases with linear interpolation between adjacent phases
- Apply anti-aliasing cutoff at min(in_rate, out_rate) / 2
- Scale filter coefficients by min(1, out_rate / in_rate) when downsampling to preserve energy
- Handle signal boundaries by clamping input indices to [0, input_len-1]
- Allocate the polyphase table internally (malloc) and free before returning
- Return the number of samples written to output
- Assert: input_len > 0, in_rate > 0, out_rate > 0, num_zero_crossings > 0, max_output_len >= MD_resample_output_len(input_len, in_rate, out_rate)

#### Scenario: Identity resampling (1:1 ratio)
- **WHEN** a signal is resampled with in_rate == out_rate
- **THEN** the output SHALL match the input within 1e-6 tolerance per sample

#### Scenario: DC signal preservation
- **WHEN** a constant signal (all 1.0) is resampled at any rate ratio
- **THEN** all output samples SHALL be within 1e-6 of 1.0

#### Scenario: Sine wave frequency preservation (44100 to 48000)
- **WHEN** a 440 Hz sine wave at 44100 Hz is resampled to 48000 Hz with num_zero_crossings=32, kaiser_beta=10.0
- **THEN** the dominant frequency in the output SHALL be 440 Hz within 1 Hz tolerance

#### Scenario: Output length matches prediction
- **WHEN** a signal is resampled for rate pairs (44100→48000, 48000→44100, 16000→8000, 8000→16000, 44100→22050)
- **THEN** the returned sample count SHALL equal `MD_resample_output_len(input_len, in_rate, out_rate)` for each pair

#### Scenario: Energy preservation
- **WHEN** white noise is resampled at any rate ratio with num_zero_crossings=32, kaiser_beta=10.0
- **THEN** the RMS of the output SHALL be within 0.5 dB of the RMS of the input

#### Scenario: Anti-aliasing suppression
- **WHEN** a 20000 Hz sine at 48000 Hz sample rate is downsampled to 16000 Hz
- **THEN** the output spectrum SHALL show no energy above -60 dB at the aliased frequency relative to the original amplitude

#### Scenario: Common rate pair conversions
- **WHEN** signals are resampled between common rate pairs (44100↔48000, 44100↔22050, 48000↔16000)
- **THEN** spectral content below the target Nyquist SHALL be preserved and content above SHALL be suppressed
