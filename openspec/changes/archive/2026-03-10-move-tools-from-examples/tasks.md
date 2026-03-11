## 1. Create tools/audio_steg directory and Makefile

- [x] 1.1 Create `tools/audio_steg/` directory
- [x] 1.2 Move `examples/audio_steg.c` to `tools/audio_steg/audio_steg.c`
- [x] 1.3 Move `examples/space_invader.png` to `tools/audio_steg/space_invader.png`
- [x] 1.4 Move `examples/minidsp_qr.png` to `tools/audio_steg/minidsp_qr.png`
- [x] 1.5 Create `tools/audio_steg/Makefile` following the mel_viz pattern (include ../../config.mk, .PHONY library delegation, link flags for fftw3, sndfile, math)

## 2. Remove audio_steg from examples build

- [x] 2.1 Remove `audio_steg` from `SNDFILE_EXAMPLES` in `examples/Makefile`
- [x] 2.2 Remove `./audio_steg` from the `plot` target in `examples/Makefile`

## 3. Update root Makefile

- [x] 3.1 Add `$(MAKE) -C tools/audio_steg` to the `tools` target
- [x] 3.2 Add `$(MAKE) -C tools/audio_steg clean` to the `clean` target

## 4. Update Doxyfile

- [x] 4.1 Add `tools/audio_steg` to `EXAMPLE_PATH` (for \snippet resolution)
- [x] 4.2 Update `HTML_EXTRA_FILES` entries: `examples/space_invader.png` → `tools/audio_steg/space_invader.png`, `examples/minidsp_qr.png` → `tools/audio_steg/minidsp_qr.png`

## 5. Update guide documentation

- [x] 5.1 Update `guides/audio-steganography.md` prose paths from `examples/audio_steg.c` to `tools/audio_steg/audio_steg.c`
- [x] 5.2 Update build commands in the guide from `make -C examples audio_steg` to `make -C tools/audio_steg`

## 6. Update ignore files

- [x] 6.1 Remove `examples/audio_steg`, `examples/audio_steg.csv`, `examples/audio_steg.html` from `.gitignore`
- [x] 6.2 Add `tools/audio_steg/audio_steg` to `.gitignore`
- [x] 6.3 Update `.dockerignore` if it has audio_steg entries (none found — no change needed)

## 7. Verify

- [x] 7.1 Run `make tools` from repo root — both mel_viz and audio_steg build
- [x] 7.2 Run `make -C examples` — builds succeed without audio_steg
- [x] 7.3 Run `make clean` from repo root — cleans audio_steg artifacts
- [x] 7.4 Grep for any remaining `examples/audio_steg` or `examples/space_invader` references (fixed `gen_audio_samples.c` path; remaining refs are only in openspec change docs)
