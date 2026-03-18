## 1. Project Setup

- [x] 1.1 Create `tools/resample/` directory
- [x] 1.2 Create `tools/resample/Makefile` following the `mel_viz`/`audio_steg` pattern (`include ../../config.mk`, `.PHONY` library delegation, link `-lminidsp -lfftw3 -lsndfile -lm`)
- [x] 1.3 Add `tools/resample/resample` to `.gitignore` and `.dockerignore`

## 2. Core Implementation

- [x] 2.1 Create `tools/resample/resample.c` with argument parsing: optional `-z`/`-b` flags, then three positional args (input file, target rate, output file)
- [x] 2.2 Implement input validation: check argument count, reject zero/negative target rate, warn on non-`.wav` output extension
- [x] 2.3 Read input audio via `FIO_read_audio()`, convert `float*` to `double*`
- [x] 2.4 Handle same-rate detection: if input rate equals target rate, copy directly and print note
- [x] 2.5 Compute output length via `MD_resample_output_len()`, allocate output buffer, call `MD_resample()`
- [x] 2.6 Convert resampled `double*` back to `float*`, write via `FIO_write_wav()`
- [x] 2.7 Print conversion summary (input/output file names, sample rates, sample counts, duration)
- [x] 2.8 Free all allocated memory and return appropriate exit code

## 3. Build Integration

- [x] 3.1 Add `$(MAKE) -C tools/resample` to root Makefile `tools` target
- [x] 3.2 Add `$(MAKE) -C tools/resample clean` to root Makefile `clean` target

## 4. Verification

- [x] 4.1 Build with `make tools` from repo root and verify binary is produced
- [x] 4.2 Test with a sample WAV file: resample to a different rate, verify output plays correctly
- [x] 4.3 Test error cases: missing args, zero rate, non-existent input file
