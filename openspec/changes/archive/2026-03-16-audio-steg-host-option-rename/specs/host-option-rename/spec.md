## ADDED Requirements

### Requirement: Usage text displays AUDIO_HOST metavar
The `usage()` output SHALL display `-i AUDIO_HOST` instead of `-i HOST` in all usage lines where the `-i` option appears.

#### Scenario: Encode usage line
- **WHEN** user runs `audio_steg -h` or triggers usage output
- **THEN** the encode usage line reads `./audio_steg --encode METHOD MSG [-i AUDIO_HOST] [-o OUT]`

#### Scenario: Encode-image usage line
- **WHEN** user runs `audio_steg -h` or triggers usage output
- **THEN** the encode-image usage line reads `./audio_steg --encode-image METHOD IMAGE [-i AUDIO_HOST] [-o OUT]`

### Requirement: File-header comment matches usage text
The file-header comment block in `audio_steg.c` SHALL use `AUDIO_HOST` (not `HOST`) in the mode descriptions for `--encode` and `--encode-image`.

#### Scenario: Header comment updated
- **WHEN** a developer reads the top-of-file comment block
- **THEN** both encode mode descriptions show `-i AUDIO_HOST.wav` instead of `-i HOST.wav`

### Requirement: Examples section includes -i usage
The "Examples:" section of the help output SHALL include at least one example demonstrating the `-i AUDIO_HOST` option.

#### Scenario: Example with -i flag present
- **WHEN** user views help output
- **THEN** there is an example line showing `--encode` with `-i host.wav` (or similar audio filename) and `-o stego.wav`
