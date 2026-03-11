## ADDED Requirements

### Requirement: audio_steg lives in tools directory
The `audio_steg` program SHALL reside at `tools/audio_steg/audio_steg.c` with its asset files (`space_invader.png`, `minidsp_qr.png`) in the same directory.

#### Scenario: Directory structure matches mel_viz pattern
- **WHEN** a developer lists the `tools/audio_steg/` directory
- **THEN** it SHALL contain `Makefile`, `audio_steg.c`, `space_invader.png`, and `minidsp_qr.png`

### Requirement: audio_steg has its own Makefile
The `tools/audio_steg/Makefile` SHALL follow the `tools/mel_viz/Makefile` pattern: include `../../config.mk`, delegate library builds to the root Makefile via `.PHONY`, and link against `libminidsp.a`, `fftw3`, `sndfile`, and `math`.

#### Scenario: Build audio_steg from its directory
- **WHEN** a developer runs `make` in `tools/audio_steg/`
- **THEN** the `audio_steg` binary SHALL be produced

#### Scenario: Build via root Makefile tools target
- **WHEN** a developer runs `make tools` from the repo root
- **THEN** `audio_steg` SHALL be built alongside `mel_viz`

### Requirement: audio_steg removed from examples
The `audio_steg` entry SHALL be removed from `examples/Makefile` (both the `SNDFILE_EXAMPLES` variable and the `plot` target). The files `examples/audio_steg.c`, `examples/space_invader.png`, and `examples/minidsp_qr.png` SHALL no longer exist.

#### Scenario: Examples directory has no audio_steg artifacts
- **WHEN** a developer runs `make -C examples` and lists the examples directory
- **THEN** no `audio_steg` binary, source, or asset files SHALL be present

### Requirement: Doxygen snippet resolution works after move
The `EXAMPLE_PATH` in `Doxyfile` SHALL include `tools/audio_steg` so that `\snippet audio_steg.c <id>` directives in `guides/audio-steganography.md` resolve correctly.

#### Scenario: Doxygen build succeeds with snippets
- **WHEN** a developer runs `make docs`
- **THEN** the audio steganography guide SHALL render with all snippet code blocks intact and no Doxygen warnings about missing snippet files

### Requirement: Documentation paths updated
All human-readable path references to `examples/audio_steg.c` in `guides/audio-steganography.md` SHALL be updated to `tools/audio_steg/audio_steg.c`. Build commands SHALL reference the new location (e.g., `make -C tools/audio_steg`).

#### Scenario: Guide references correct paths
- **WHEN** a developer reads `guides/audio-steganography.md`
- **THEN** all file paths and build commands SHALL reference `tools/audio_steg/`, not `examples/`

### Requirement: HTML_EXTRA_FILES updated for moved assets
The `Doxyfile` `HTML_EXTRA_FILES` entries for `examples/space_invader.png` and `examples/minidsp_qr.png` SHALL be updated to `tools/audio_steg/space_invader.png` and `tools/audio_steg/minidsp_qr.png`.

#### Scenario: Doxygen copies PNG assets correctly
- **WHEN** Doxygen generates HTML documentation
- **THEN** `space_invader.png` and `minidsp_qr.png` SHALL be copied to the output directory from their new paths

### Requirement: Git ignore and Docker ignore updated
`.gitignore` entries for `examples/audio_steg`, `examples/audio_steg.csv`, and `examples/audio_steg.html` SHALL be removed. Appropriate entries for `tools/audio_steg/audio_steg` (the binary) SHALL be added.

#### Scenario: Build artifacts are ignored
- **WHEN** a developer builds `audio_steg` in `tools/audio_steg/`
- **THEN** the binary and any `.dSYM` directory SHALL be ignored by git

### Requirement: Root Makefile clean target updated
The root Makefile `clean` target SHALL include `$(MAKE) -C tools/audio_steg clean`.

#### Scenario: Clean removes audio_steg artifacts
- **WHEN** a developer runs `make clean` from the repo root
- **THEN** the `tools/audio_steg/audio_steg` binary and `.dSYM` directory SHALL be removed
