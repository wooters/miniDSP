## Context

The `audio_steg` tool's `-i` option uses the metavar `HOST` in usage text. This is unclear — it could be confused with a network host. The actual argument is an audio WAV file that serves as the carrier signal. Renaming to `AUDIO_HOST` and adding an example will make the CLI self-documenting.

## Goals / Non-Goals

**Goals:**
- Rename `-i HOST` to `-i AUDIO_HOST` in all usage/help text
- Add an example line demonstrating `-i AUDIO_HOST` usage
- Update the file-header comment block to match

**Non-Goals:**
- Changing the actual `-i` flag name (it stays `-i`)
- Modifying any runtime behavior or argument parsing logic
- Renaming variables in the C code (e.g., `infile` stays `infile`)

## Decisions

1. **Metavar name `AUDIO_HOST`** — Descriptive without being verbose. Alternatives considered: `HOST_WAV` (too format-specific), `CARRIER` (correct DSP term but less intuitive for casual users), `INPUT` (conflicts with the general meaning of `-i`).

2. **Example placement** — Add the new example after the first `--encode` example (which omits `-i`), so users see the basic form first, then the variant with a custom host file.

3. **Example wording** — Use `--encode lsb "secret message" -i host.wav -o stego.wav` to show both `-i` and `-o` together, since that's the most common real-world usage.

## Risks / Trade-offs

- None. This is a cosmetic change to help text only.
