## ADDED Requirements

### Requirement: Windowed-sinc FIR lowpass filter design
The library SHALL provide `MD_design_lowpass_fir(double *coeffs, unsigned num_taps, double cutoff_freq, double sample_rate, double kaiser_beta)` generating a Kaiser-windowed sinc lowpass filter.

The function SHALL:
- Compute normalized cutoff: fc = cutoff_freq / sample_rate
- For each tap i: coeffs[i] = 2 * fc * sinc(2 * fc * (i - (num_taps-1)/2.0)) * kaiser[i]
- Normalize coefficients to sum to 1.0 (unity DC gain)
- Assert: num_taps > 0, sample_rate > 0, cutoff_freq > 0, cutoff_freq < sample_rate / 2

#### Scenario: DC gain is unity
- **WHEN** `MD_design_lowpass_fir` generates a filter with any valid parameters
- **THEN** the sum of all coefficients SHALL be within 1e-12 of `1.0`

#### Scenario: Linear phase symmetry
- **WHEN** a lowpass filter is designed with any valid parameters
- **THEN** `coeffs[i]` SHALL equal `coeffs[num_taps-1-i]` for all valid i (within floating-point tolerance)

#### Scenario: Passband response
- **WHEN** a lowpass filter (cutoff = 4000 Hz, sample_rate = 48000 Hz, 65 taps, β = 10.0) is convolved with a 1000 Hz sine wave at 48000 Hz
- **THEN** the output amplitude SHALL be within 0.1 dB of the input amplitude (passband flatness)

#### Scenario: Stopband rejection
- **WHEN** a lowpass filter (cutoff = 4000 Hz, sample_rate = 48000 Hz, 65 taps, β = 10.0) is convolved with a 10000 Hz sine wave at 48000 Hz
- **THEN** the output amplitude SHALL be at least 60 dB below the input amplitude

#### Scenario: Assertions on invalid parameters
- **WHEN** `MD_design_lowpass_fir` is called with cutoff_freq >= sample_rate / 2
- **THEN** the function SHALL trigger an assertion failure
