## ADDED Requirements

### Requirement: CLI accepts input file, target rate, and output file
The tool SHALL accept three positional arguments: input file path, target sample rate (integer Hz), and output file path.

#### Scenario: Valid invocation
- **WHEN** user runs `resample input.wav 48000 output.wav`
- **THEN** the tool reads `input.wav`, resamples to 48000 Hz, and writes `output.wav`

#### Scenario: Missing arguments
- **WHEN** user runs `resample` with fewer than 3 positional arguments
- **THEN** the tool prints a usage message to stderr and exits with non-zero status

### Requirement: Tool reads mono audio files only
The tool SHALL read single-channel (mono) audio files using `FIO_read_audio()`. If the input file has more than one channel, the tool SHALL print an error message and exit with non-zero status.

#### Scenario: Mono WAV input
- **WHEN** user provides a mono WAV file as input
- **THEN** the tool reads it successfully

#### Scenario: Multi-channel input rejected
- **WHEN** user provides a stereo or multi-channel audio file
- **THEN** the tool prints "Error: only mono (single-channel) files are supported" to stderr and exits with non-zero status

### Requirement: Tool writes WAV output only
The tool SHALL write output exclusively in WAV format (IEEE float) using `FIO_write_wav()`. The usage message SHALL state that only `.wav` output is supported.

#### Scenario: WAV output written
- **WHEN** resampling succeeds
- **THEN** the tool writes an IEEE float WAV file to the specified output path

#### Scenario: Non-WAV extension warning
- **WHEN** the output filename does not end in `.wav`
- **THEN** the tool prints a warning to stderr that the output will be WAV format regardless of extension, and proceeds with writing

### Requirement: Resampling uses sensible defaults
The tool SHALL use 32 zero-crossings and kaiser beta 10.0 as default resampler parameters, matching the `MD_resample()` documentation examples.

#### Scenario: Default quality parameters
- **WHEN** user runs `resample input.wav 48000 output.wav` without `-z` or `-b` flags
- **THEN** the tool calls `MD_resample()` with `num_zero_crossings=32` and `kaiser_beta=10.0`

### Requirement: Optional quality parameter overrides
The tool SHALL accept optional `-z <zero_crossings>` and `-b <kaiser_beta>` flags to override the default resampler quality parameters.

#### Scenario: Custom zero-crossings
- **WHEN** user runs `resample -z 64 input.wav 48000 output.wav`
- **THEN** the tool calls `MD_resample()` with `num_zero_crossings=64`

#### Scenario: Custom kaiser beta
- **WHEN** user runs `resample -b 14.0 input.wav 48000 output.wav`
- **THEN** the tool calls `MD_resample()` with `kaiser_beta=14.0`

### Requirement: Tool prints conversion summary
The tool SHALL print a summary to stdout showing input file name, input sample rate, input duration, output sample rate, output sample count, and output file name.

#### Scenario: Summary output
- **WHEN** resampling completes successfully
- **THEN** the tool prints a summary like:
  ```
  Input:  input.wav (44100 Hz, 132300 samples, 3.00 s)
  Output: output.wav (48000 Hz, 144000 samples)
  ```

### Requirement: Same-rate detection
The tool SHALL detect when the input sample rate equals the target sample rate, copy the audio directly without resampling, and print a note that no resampling was needed.

#### Scenario: Input rate equals target rate
- **WHEN** user runs `resample input.wav 44100 output.wav` and `input.wav` is already 44100 Hz
- **THEN** the tool writes an identical copy and prints "Note: input is already at 44100 Hz, copying without resampling"

### Requirement: Invalid target rate rejected
The tool SHALL reject target sample rates that are zero or negative.

#### Scenario: Zero target rate
- **WHEN** user runs `resample input.wav 0 output.wav`
- **THEN** the tool prints an error to stderr and exits with non-zero status

### Requirement: Build integration
The tool SHALL be built by `make tools` from the repo root and cleaned by `make clean`.

#### Scenario: Build via make tools
- **WHEN** user runs `make tools` from the repo root
- **THEN** `tools/resample/resample` binary is built

#### Scenario: Clean removes binary
- **WHEN** user runs `make clean` from the repo root
- **THEN** `tools/resample/resample` and its `.dSYM` directory are removed
