## Why

miniDSP already provides `MD_resample()` for polyphase sinc resampling, but there is no command-line tool that exposes this capability. Users who want to convert a WAV file from one sample rate to another must write their own C program. A standalone CLI tool in `tools/` makes sample rate conversion a one-liner, following the same pattern as `mel_viz` and `audio_steg`.

## What Changes

- Add a new CLI tool `tools/resample/resample` that reads a mono audio file, resamples it to a target rate, and writes a WAV output file.
- Accept three required command-line arguments: input file path, target sample rate, and output file path.
- Validate inputs and print clear error messages (mono-only constraint, WAV-only output).
- Integrate into the root Makefile `tools` target and clean rules.
- Update `.gitignore` and `.dockerignore` for the new binary.

## Capabilities

### New Capabilities
- `resample-cli`: Command-line interface for sample rate conversion — argument parsing, validation, invocation of `MD_resample()`, and file I/O via `FIO_read_audio()` / `FIO_write_wav()`.

### Modified Capabilities
(none)

## Impact

- **New files**: `tools/resample/resample.c`, `tools/resample/Makefile`
- **Modified files**: root `Makefile` (tools/clean targets), `.gitignore`, `.dockerignore`
- **Dependencies**: No new external dependencies — uses existing `libminidsp.a`, `libsndfile`, `libfftw3`
- **API**: No library API changes
