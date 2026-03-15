## ADDED Requirements

### Requirement: MD_adjust_dblevel handles zero-energy input
`MD_adjust_dblevel` SHALL produce correct output when the input signal has zero energy (all-zero or near-zero samples), instead of producing NaN/Inf.

#### Scenario: Silent input produces silent output
- **WHEN** `MD_adjust_dblevel` is called with an all-zero input signal
- **THEN** the output signal SHALL be all zeros (silence copied unchanged)

#### Scenario: Near-zero input produces finite output
- **WHEN** `MD_adjust_dblevel` is called with an extremely quiet signal (energy approaching zero)
- **THEN** all output samples SHALL be finite (no NaN or Inf values)

### Requirement: BiQuad_new rejects degenerate frequencies
`BiQuad_new` SHALL return NULL when the requested filter frequency would cause a division by zero in coefficient calculation.

#### Scenario: Frequency of zero returns NULL
- **WHEN** `BiQuad_new` is called with `freq = 0.0`
- **THEN** it SHALL return NULL

#### Scenario: Frequency at Nyquist returns NULL
- **WHEN** `BiQuad_new` is called with `freq = srate / 2.0`
- **THEN** it SHALL return NULL

#### Scenario: Frequencies just inside valid range succeed
- **WHEN** `BiQuad_new` is called with `freq = 1.0` (just above 0) and `freq = srate/2 - 1` (just below Nyquist)
- **THEN** it SHALL return a valid non-NULL biquad pointer
