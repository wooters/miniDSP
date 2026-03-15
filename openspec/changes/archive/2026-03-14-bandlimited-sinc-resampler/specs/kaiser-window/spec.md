## ADDED Requirements

### Requirement: Kaiser window generator
The library SHALL provide `MD_Gen_Kaiser_Win(double *out, unsigned n, double beta)` generating a Kaiser window using the formula w[i] = I₀(β √(1 - ((2i/(n-1)) - 1)²)) / I₀(β).

The function SHALL assert `n > 0` and follow the existing `MD_Gen_*_Win` family conventions.

#### Scenario: Degenerate single-sample window
- **WHEN** `MD_Gen_Kaiser_Win(out, 1, 10.0)` is called
- **THEN** `out[0]` SHALL be `1.0`

#### Scenario: Window symmetry
- **WHEN** a Kaiser window of length n > 1 is generated with any β > 0
- **THEN** `out[i]` SHALL equal `out[n-1-i]` for all valid i (within floating-point tolerance)

#### Scenario: Tapered ends with peak at center
- **WHEN** a Kaiser window of length n ≥ 3 is generated with β > 0
- **THEN** `out[0]` SHALL be greater than 0, and `out[0]` SHALL be less than `out[n/2]` (the center value)

#### Scenario: β controls sidelobe attenuation
- **WHEN** Kaiser windows are generated with β = 5.0 and β = 14.0 (same length)
- **THEN** the β = 14.0 window SHALL have a smaller ratio of `out[0] / out[n/2]` (more tapered ends relative to center)

#### Scenario: Known β values produce expected stopband attenuation
- **WHEN** a Kaiser window with β ≈ 10.0 is used in a windowed-sinc lowpass filter
- **THEN** the stopband attenuation SHALL be approximately 100 dB
