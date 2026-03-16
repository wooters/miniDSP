## ADDED Requirements

### Requirement: Brickwall lowpass function signature
The library SHALL provide a public function:
```c
void MD_lowpass_brickwall(double *signal, unsigned len, double cutoff_hz, double sample_rate);
```
The function SHALL modify `signal` in-place, zeroing all frequency content above `cutoff_hz`.

#### Scenario: Function exists and compiles
- **WHEN** a program includes `"minidsp.h"` and calls `MD_lowpass_brickwall(buf, 1024, 8000.0, 48000.0)`
- **THEN** it SHALL compile and link without errors

### Requirement: Brickwall lowpass filtering behavior
`MD_lowpass_brickwall()` SHALL:
1. Compute an r2c FFT of the input signal
2. Zero all complex bins whose center frequency exceeds `cutoff_hz` (i.e., bins with index > `floor(cutoff_hz * len / sample_rate)`)
3. Compute a c2r IFFT of the modified spectrum
4. Normalize the result by dividing each sample by `len` (to compensate for FFTW's unnormalized transform)
5. Write the result back to the `signal` buffer

The function SHALL preserve signal content at and below `cutoff_hz` and completely eliminate content above `cutoff_hz`.

#### Scenario: Pure tone below cutoff is preserved
- **WHEN** a 1000 Hz sine wave at amplitude 1.0 is generated at 48000 Hz sample rate
- **AND** `MD_lowpass_brickwall()` is applied with `cutoff_hz = 8000.0`
- **THEN** the output signal energy SHALL be within 0.1% of the input signal energy

#### Scenario: Pure tone above cutoff is eliminated
- **WHEN** a 20000 Hz sine wave at amplitude 1.0 is generated at 48000 Hz sample rate
- **AND** `MD_lowpass_brickwall()` is applied with `cutoff_hz = 8000.0`
- **THEN** the output signal energy SHALL be less than 1e-20 (effectively zero)

#### Scenario: Mixed signal retains only low frequencies
- **WHEN** a signal containing both a 1000 Hz tone and a 20000 Hz tone is created at 48000 Hz
- **AND** `MD_lowpass_brickwall()` is applied with `cutoff_hz = 8000.0`
- **THEN** the output SHALL contain only the 1000 Hz tone
- **AND** the energy of the 20000 Hz component SHALL be less than 1e-20

### Requirement: Brickwall lowpass preconditions
`MD_lowpass_brickwall()` SHALL assert the following preconditions:
- `signal` is not NULL
- `len` is greater than 0
- `cutoff_hz` is greater than 0
- `cutoff_hz` is less than `sample_rate / 2`
- `sample_rate` is greater than 0

#### Scenario: NULL signal assertion
- **WHEN** `MD_lowpass_brickwall(NULL, 1024, 8000.0, 48000.0)` is called
- **THEN** the program SHALL abort via assertion failure

### Requirement: Brickwall lowpass Doxygen documentation
`include/minidsp.h` SHALL contain a Doxygen doc-comment for `MD_lowpass_brickwall()` that includes:
- A `@brief` description
- The brickwall filter formula using `\f[...\f]` display math
- `@param` documentation for all four parameters
- A `@code` usage example

The feature list at the top of `minidsp.h` SHALL include brickwall lowpass filtering.

#### Scenario: Feature list updated
- **WHEN** a user reads the header comment in `minidsp.h`
- **THEN** the bulleted feature list SHALL mention FFT-based brickwall lowpass filtering

### Requirement: Brickwall lowpass tests
`tests/test_minidsp_spectrum.c` SHALL include tests for `MD_lowpass_brickwall()` covering:
- A tone below the cutoff is preserved (energy unchanged)
- A tone above the cutoff is eliminated (energy ~ 0)
- A mixed-frequency signal retains only components below the cutoff

Tests SHALL follow the project convention: `static int test_xxx(void)` returning 1=pass, 0=fail, registered with `RUN_TEST()`.

#### Scenario: Tests pass
- **WHEN** `make test` is run
- **THEN** all `MD_lowpass_brickwall` tests SHALL pass
