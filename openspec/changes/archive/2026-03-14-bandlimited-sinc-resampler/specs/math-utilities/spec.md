## ADDED Requirements

### Requirement: Bessel I₀ function
The library SHALL provide `MD_bessel_i0(double x)` returning the zeroth-order modified Bessel function of the first kind, computed via the power series I₀(x) = Σ [(x/2)^k / k!]² until the term is below 1e-15 relative tolerance.

The function SHALL be a pure function with no allocations or side effects.

#### Scenario: I₀(0) equals 1
- **WHEN** `MD_bessel_i0(0.0)` is called
- **THEN** the return value SHALL be exactly `1.0`

#### Scenario: Known value at x=1
- **WHEN** `MD_bessel_i0(1.0)` is called
- **THEN** the return value SHALL be within 1e-10 of `1.2660658777520084`

#### Scenario: Known value at x=5
- **WHEN** `MD_bessel_i0(5.0)` is called
- **THEN** the return value SHALL be within 1e-6 of `27.239871823604442`

#### Scenario: Monotonic increase for positive x
- **WHEN** `MD_bessel_i0` is evaluated at x = 0, 1, 2, 3, 4, 5
- **THEN** each successive return value SHALL be strictly greater than the previous

### Requirement: Normalized sinc function
The library SHALL provide `MD_sinc(double x)` returning sin(πx)/(πx), with a near-zero threshold of |x| < 1e-12 returning 1.0.

#### Scenario: sinc(0) equals 1
- **WHEN** `MD_sinc(0.0)` is called
- **THEN** the return value SHALL be exactly `1.0`

#### Scenario: sinc at integer values equals 0
- **WHEN** `MD_sinc(n)` is called for integer n in {-3, -2, -1, 1, 2, 3}
- **THEN** the return value SHALL be within 1e-15 of `0.0`

#### Scenario: sinc(0.5) equals 2/π
- **WHEN** `MD_sinc(0.5)` is called
- **THEN** the return value SHALL be within 1e-12 of `2.0 / M_PI`

#### Scenario: Near-zero threshold
- **WHEN** `MD_sinc(1e-13)` is called
- **THEN** the return value SHALL be `1.0` (handled by the near-zero branch)
