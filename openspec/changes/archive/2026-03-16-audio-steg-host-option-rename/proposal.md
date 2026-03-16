## Why

The `-i HOST` option in `audio_steg` is ambiguous — "HOST" could be mistaken for a network hostname. Renaming it to `-i AUDIO_HOST` makes it clear the argument is an audio file used as the carrier signal. Additionally, none of the existing help examples demonstrate the `-i` flag, so users may not know how to use a custom host file.

## What Changes

- Rename the `-i HOST` metavar to `-i AUDIO_HOST` in all usage lines
- Add an example to the "Examples:" section showing `-i AUDIO_HOST` in use with `--encode`

## Capabilities

### New Capabilities
- `host-option-rename`: Rename `-i HOST` metavar to `-i AUDIO_HOST` in usage text and add a usage example demonstrating the flag

### Modified Capabilities

## Impact

- `tools/audio_steg/audio_steg.c` — `usage()` function (usage lines + examples section)
- The file-header comment block (lines 5–26) that documents the CLI modes
- No behavioral or API changes; cosmetic/documentation only
