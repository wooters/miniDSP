## Context

`examples/audio_steg.c` is a 641-line, 6-mode CLI steganography program that embeds assets (`space_invader.png`, `minidsp_qr.png`). It sits alongside simple API demos like `sine_wave.c` (181 lines) and `impulse.c` (195 lines). The `tools/` directory was created for exactly this kind of substantial program — `tools/mel_viz/` (460 lines + web UI) is the established pattern.

The Doxygen guide `guides/audio-steganography.md` references `examples/audio_steg.c` via `\snippet` directives, and `Doxyfile` lists `examples/space_invader.png` in `HTML_EXTRA_FILES` and uses `EXAMPLE_PATH = examples` for snippet resolution.

## Goals / Non-Goals

**Goals:**
- Move `audio_steg.c`, `space_invader.png`, and `minidsp_qr.png` to `tools/audio_steg/`
- Give `audio_steg` its own Makefile following the `mel_viz` pattern
- Update all build system, documentation, and ignore-file references
- Keep `\snippet` directives in `guides/audio-steganography.md` working

**Non-Goals:**
- Refactoring `audio_steg.c` code itself — pure move, no logic changes
- Moving any other examples (this change is scoped to `audio_steg` only)
- Changing the library API or test suite

## Decisions

**1. Directory structure follows mel_viz pattern**

```
tools/audio_steg/
├── Makefile          (include ../../config.mk)
├── audio_steg.c
├── space_invader.png
└── minidsp_qr.png
```

*Rationale*: Consistency with the existing tool. Same include paths, same library delegation pattern.

**2. Add `tools/audio_steg` to Doxyfile `EXAMPLE_PATH`**

Currently `EXAMPLE_PATH = examples`. Snippet directives use bare filenames (`\snippet audio_steg.c self-test`), which resolve by searching `EXAMPLE_PATH` directories. Adding `tools/audio_steg` to `EXAMPLE_PATH` makes snippets resolve without changing any guide markdown.

*Alternative considered*: Changing `\snippet` paths to absolute — rejected because it's fragile and non-standard.

**3. Update `HTML_EXTRA_FILES` paths for PNG assets**

Change `examples/space_invader.png` and `examples/minidsp_qr.png` to `tools/audio_steg/space_invader.png` and `tools/audio_steg/minidsp_qr.png`.

**4. Update guide prose paths but not snippet directives**

The guide has human-readable text like `examples/audio_steg.c` and `make -C examples audio_steg`. These need updating. The `\snippet audio_steg.c id` directives do NOT need updating (Doxygen resolves by filename search in `EXAMPLE_PATH`).

## Risks / Trade-offs

- **[Broken snippet resolution]** → Mitigated by adding `tools/audio_steg` to `EXAMPLE_PATH` and verifying with `make docs`
- **[Stale references]** → Grep for `audio_steg` and `space_invader` across all files to catch every reference
- **[Build regression]** → Verify `make tools` and `make clean` both work after the move
