## ADDED Requirements

### Requirement: Platform compatibility note
The README SHALL include a note near the top stating the library compiles and runs on Ubuntu and macOS, that Windows is untested, and that pull requests for Windows support are welcome.

#### Scenario: Reader sees platform info
- **WHEN** a user opens the README
- **THEN** they see supported platforms (Ubuntu, macOS) and a note about Windows before any code or build instructions

### Requirement: audio_steg tool listing
The "Tools" section SHALL include an entry for `audio_steg` describing it as an audio steganography tool that hides and recovers secret messages or binary data in WAV files.

#### Scenario: Tools section lists audio_steg
- **WHEN** a user reads the "Tools" section
- **THEN** they see entries for both `mel_viz` and `audio_steg` with descriptions and build/run instructions

### Requirement: Build section before Use section
The "Build and Test" section (including all subsections) SHALL appear before the "Use in your project" section in the README.

#### Scenario: Section ordering
- **WHEN** a user reads the README top to bottom
- **THEN** they encounter "Build and Test" before "Use in your project"

### Requirement: Accurate Apple container dependency info
The Dependencies table SHALL state that the Apple `container` CLI must be installed from its GitHub repository, not that it is "built-in" on macOS 26+.

#### Scenario: Container row accuracy
- **WHEN** a user reads the Dependencies table
- **THEN** the Apple container row says to install from GitHub (with link), not "built-in"
