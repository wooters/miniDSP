## ADDED Requirements

### Requirement: Spectext method constant
The library SHALL define `MD_STEG_SPECTEXT` with integer value `2`.

### Requirement: Spectext capacity calculation
`MD_steg_capacity(signal_len, sample_rate, MD_STEG_SPECTEXT)` SHALL return the maximum number of message bytes that can be encoded.

The capacity SHALL be the minimum of:
- The LSB capacity of the (possibly resampled) signal: `(out_len - 32) / 8`
- The visual capacity based on host duration: `floor(duration_sec / 0.24)` characters

where `out_len` is `MD_resample_output_len(signal_len, sample_rate, 48000.0)` if `sample_rate < 48000`, otherwise `signal_len`, and `duration_sec = signal_len / sample_rate`.

#### Scenario: Capacity of a 3-second 44.1 kHz signal
- **WHEN** `MD_steg_capacity(132300, 44100.0, MD_STEG_SPECTEXT)` is called
- **THEN** the return value SHALL be `12` (limited by visual capacity: floor(3.0 / 0.24) = 12)

#### Scenario: Capacity of a 30-second 48 kHz signal
- **WHEN** `MD_steg_capacity(1440000, 48000.0, MD_STEG_SPECTEXT)` is called
- **THEN** the return value SHALL be `125` (limited by visual capacity: floor(30.0 / 0.24) = 125)

### Requirement: Spectext encode
`MD_steg_encode(host, output, signal_len, sample_rate, message, MD_STEG_SPECTEXT)` SHALL:
1. Upsample the host signal to 48 kHz using `MD_resample()` if `sample_rate < 48000.0`
2. LSB-encode the full message into the (resampled) host using the existing LSB framing (32-bit header + payload)
3. Compute the number of visually displayable characters: `vis_chars = min(strlen(message), floor(duration_sec / 0.24))` where `duration_sec = signal_len / sample_rate`
4. Generate spectrogram text art for the first `vis_chars` characters using `MD_spectrogram_text()` with `freq_lo=18000.0`, `freq_hi=24000.0`, fixed 30 ms column width, and `sample_rate=48000.0`
5. Scale the spectrogram text signal to 0.02 peak amplitude (multiply by `0.02 / 0.9`)
6. Additively mix the scaled spectrogram text into the LSB-encoded output
7. Write the result to the `output` buffer (which may be longer than `signal_len` if upsampling occurred)
8. Return the number of message bytes encoded

The function SHALL assert that the message length does not exceed the LSB capacity of the output signal.

#### Scenario: Round-trip encode and decode
- **WHEN** a message "HELLO" is encoded with `MD_STEG_SPECTEXT` into a 3-second 48 kHz host signal
- **AND** decoded with `MD_steg_decode(..., MD_STEG_SPECTEXT)`
- **THEN** the decoded message SHALL be "HELLO"

#### Scenario: Upsampling from 44.1 kHz
- **WHEN** a message is encoded with `MD_STEG_SPECTEXT` into a host signal at 44100 Hz
- **THEN** the output signal length SHALL be `MD_resample_output_len(signal_len, 44100.0, 48000.0)`

#### Scenario: No upsampling at 48 kHz
- **WHEN** a message is encoded with `MD_STEG_SPECTEXT` into a host signal already at 48000 Hz
- **THEN** the output signal length SHALL equal `signal_len` (no resampling performed)

#### Scenario: Visual truncation for long messages
- **WHEN** a 20-character message is encoded into a 3-second host (visual capacity = 12 chars)
- **THEN** the full 20-character message SHALL be recoverable via `MD_steg_decode()`
- **AND** the spectrogram text art SHALL render only the first 12 characters

#### Scenario: Ultrasonic energy present after encoding
- **WHEN** a message is encoded with `MD_STEG_SPECTEXT`
- **THEN** the RMS energy in the 18-24 kHz band of the output SHALL be measurably above the noise floor

### Requirement: Spectext decode
`MD_steg_decode(stego, signal_len, sample_rate, message_out, max_msg_len, MD_STEG_SPECTEXT)` SHALL delegate to the existing LSB decode path and return the number of message bytes decoded.

#### Scenario: Decode matches LSB decode
- **WHEN** a spectext-encoded signal is decoded with `MD_STEG_SPECTEXT`
- **AND** the same signal is decoded with `MD_STEG_LSB`
- **THEN** both SHALL return the same message

### Requirement: Spectext binary encode/decode
`MD_steg_encode_bytes()` and `MD_steg_decode_bytes()` SHALL support `MD_STEG_SPECTEXT` with the same hybrid behavior: LSB data channel for the binary payload, spectrogram text art rendering the hex or truncated representation of the data in the ultrasonic band.

For binary payloads, the spectrogram text art SHALL render the string `"[BIN <N>B]"` where `<N>` is the payload size in bytes, since raw binary data has no meaningful text representation.

#### Scenario: Binary round-trip
- **WHEN** 100 bytes of arbitrary binary data are encoded with `MD_STEG_SPECTEXT`
- **AND** decoded with `MD_steg_decode_bytes(..., MD_STEG_SPECTEXT)`
- **THEN** the decoded data SHALL match the input exactly

### Requirement: Spectext detection
`MD_steg_detect()` SHALL return `MD_STEG_SPECTEXT` when a valid LSB framing header is present AND the signal contains significant energy in the 18-24 kHz band.

The detection order SHALL be:
1. Probe for BFSK carriers (existing)
2. Probe for LSB header (existing)
3. If LSB header found, check for ultrasonic energy in 18-24 kHz
4. If ultrasonic energy exceeds threshold, return `MD_STEG_SPECTEXT`; otherwise return `MD_STEG_LSB`

#### Scenario: Detect spectext-encoded signal
- **WHEN** `MD_steg_detect()` is called on a spectext-encoded signal
- **THEN** it SHALL return `MD_STEG_SPECTEXT`

#### Scenario: Detect plain LSB (no ultrasonic energy)
- **WHEN** `MD_steg_detect()` is called on an LSB-encoded signal (no spectrogram text)
- **THEN** it SHALL return `MD_STEG_LSB` (not `MD_STEG_SPECTEXT`)

### Requirement: Output buffer sizing
The caller SHALL provide an output buffer large enough for the 48 kHz output. For convenience, `MD_steg_capacity()` documentation SHALL note that the output buffer must accommodate `MD_resample_output_len(signal_len, sample_rate, 48000.0)` samples when `sample_rate < 48000`.

## MODIFIED Requirements

### Requirement: CLI method parsing
The `audio_steg` CLI tool SHALL accept `"spectext"` as a valid method name mapping to `MD_STEG_SPECTEXT`, in addition to existing `"lsb"` and `"freq"`.

### Requirement: CLI self-test
The `audio_steg` self-test mode (no arguments) SHALL include a spectext encode/decode round-trip test alongside the existing LSB and BFSK tests.

### Requirement: Steg guide documentation
The `guides/audio-steganography.md` tutorial SHALL be expanded with a new section covering:
- Method overview and hybrid architecture (LSB data + spectrogram art visual)
- Encode pipeline diagram
- Fixed column width and capacity formula with "Reading the formula in C" snippet
- Frequency mapping and amplitude scaling explanation
- Automatic upsampling behavior
- Listening comparison (HTML5 audio: host vs. spectext-encoded)
- Spectrogram visualization (Plotly iframe showing the embedded text in 18-24 kHz band)
- Visual truncation behavior
- Updated method comparison table (adding spectext column)
- Updated API reference with spectext-specific notes

### Requirement: Multimedia assets
The following new assets SHALL be generated and added to `Doxyfile` `HTML_EXTRA_FILES`:
- `guides/audio/steg_spectext.wav` — host signal after spectext encoding (48 kHz)
- `guides/plots/steg_spectext_spectrogram.html` — Plotly spectrogram showing embedded text in the 18-24 kHz band
