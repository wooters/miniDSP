## Context

miniDSP has two steganography methods (LSB, BFSK) in `src/minidsp_steg.c`, a spectrogram text synthesizer in `src/minidsp_spectext.c`, and a polyphase sinc resampler in `src/minidsp_resample.c`. This change composes all three into a hybrid steganography method that embeds the message via LSB (for machine decoding) and simultaneously renders it as spectrogram text art in the 18-24 kHz ultrasonic band (for visual verification).

Existing infrastructure:
- `MD_spectrogram_text()` — 5x7 bitmap font, phase-continuous sines, 3 ms crossfade, normalizes to 0.9 peak
- `MD_resample()` / `MD_resample_output_len()` — polyphase sinc, arbitrary ratios, anti-aliasing built in
- `lsb_encode()` / `lsb_decode()` — internal LSB functions in `minidsp_steg.c`, already proven
- 32-bit framing header with payload type flag — shared across all methods

## Goals / Non-Goals

**Goals:**
- Machine-readable round-trip: encode a message, decode it back identically
- Human-readable visual channel: message visible in any spectrogram viewer
- Automatic upsampling: accept any input sample rate, output 48 kHz
- Fixed-width characters: consistent visual appearance regardless of message length
- Fit cleanly into the existing `MD_steg_*` API (no new function signatures)
- Visual truncation: if message exceeds visual capacity, show first N chars; LSB carries full payload

**Non-Goals:**
- OCR-based decode from the spectrogram (decode uses LSB, not image recognition)
- Streaming / real-time encoding
- Multi-channel support (mono only, consistent with existing methods)
- Configurable font size or frequency band (hardcoded for simplicity)

## Decisions

### 1. Hybrid LSB + spectrogram art in a single method enum value

The new method `MD_STEG_SPECTEXT` (value `2`) combines LSB data encoding with spectrogram text overlay. Decode delegates to `lsb_decode()`. Detection distinguishes from plain LSB by checking for energy in the 18-24 kHz band.

**Why:** This avoids breaking the existing encode/decode API symmetry. Users get a visual bonus without losing machine-readable round-trip. The LSB channel is proven reliable; there's no reason to invent a new data channel.

### 2. Fixed 30 ms column width

Each bitmap column occupies exactly 30 ms of audio, giving ~3 spectrogram pixels per font column with typical STFT settings (N=2048, hop=512 at 48 kHz). Character width is 8 columns (5 data + 3 spacing) = 240 ms per character.

**Why:** Fixed width ensures consistent visual appearance. The text occupies only as much time as it needs — a short message leaves the rest of the signal untouched. 30 ms gives comfortable legibility without requiring specific STFT parameters to read.

**Capacity at 30 ms/col:**

| Host duration | Max visible chars |
|---------------|-------------------|
| 3 sec         | 12                |
| 10 sec        | 41                |
| 30 sec        | 125               |
| 60 sec        | 250               |

### 3. Always output 48 kHz

The spectrogram art uses the 18-24 kHz band, requiring Nyquist >= 24 kHz (sample rate >= 48 kHz). If the input is below 48 kHz, the encode function upsamples it using `MD_resample()`. The output WAV is always 48 kHz.

**Why:** Simplifies the pipeline — no conditional frequency band selection. 48 kHz is a standard audio rate supported by all modern DAWs and players. The resampler's anti-aliasing filter ensures no artifacts from upsampling.

**Alternative considered:** Require 48 kHz input and reject lower rates. Rejected because it forces the user to manually resample, adding friction.

### 4. Spectrogram text amplitude: 0.02

The spectrogram text signal is scaled to 0.02 peak amplitude (~-34 dB relative to full scale) before mixing with the host. This matches the BFSK tone amplitude.

**Why:** 0.02 is loud enough to be clearly visible in a spectrogram (well above typical noise floors of -60 to -80 dB) but quiet enough to be completely imperceptible acoustically. At 18-24 kHz, most adults cannot hear these frequencies at all, and 0.02 amplitude provides additional safety margin for younger listeners.

**Scaling approach:** `MD_spectrogram_text()` normalizes to 0.9 peak. After generation, multiply by `0.02 / 0.9 ≈ 0.0222` to reach target amplitude. This avoids modifying the spectrogram text API.

### 5. Visual truncation for long messages

If the message is longer than the visual capacity (host_duration / 0.24 characters), the spectrogram art shows the first N characters that fit. The LSB channel always carries the complete message.

**Why:** The LSB channel has vastly more capacity (~5.5 KB/sec at 48 kHz vs. ~4 chars/sec for visual). Rejecting long messages would be unnecessarily restrictive. The visual channel is a label/preview, not the primary data transport.

### 6. Detection: LSB header + ultrasonic energy probe

`MD_steg_detect()` already checks for LSB framing headers. To distinguish `MD_STEG_SPECTEXT` from plain `MD_STEG_LSB`, compute the RMS energy in the 18-24 kHz band (via a short STFT or bandpass analysis). If LSB header is valid AND ultrasonic energy exceeds a threshold, report `MD_STEG_SPECTEXT`.

**Why:** Simple, no false positives on natural audio (which has negligible energy above 18 kHz after recording/compression). The probe order is: check BFSK first (existing), then check LSB header, then check ultrasonic energy to disambiguate LSB vs. spectext.

### 7. Encode pipeline

```
spectext_encode(host, output, signal_len, sample_rate, message):

  1. Upsample host to 48 kHz if sample_rate < 48000
     - out_len = MD_resample_output_len(signal_len, sample_rate, 48000)
     - MD_resample(host, signal_len, resampled, out_len, sample_rate, 48000, 32, 10.0)

  2. LSB-encode message into resampled host
     - lsb_encode(resampled, lsb_out, out_len, message)

  3. Generate spectrogram text art
     - vis_chars = min(strlen(message), floor(host_duration / 0.24))
     - vis_text = first vis_chars characters of message
     - vis_duration = vis_chars * 8 * 0.030  (30 ms per column)
     - MD_spectrogram_text(specbuf, out_len, vis_text, 18000, 24000, vis_duration, 48000)

  4. Scale spectrogram text to steganographic amplitude
     - specbuf[i] *= 0.02 / 0.9

  5. Mix: output[i] = lsb_out[i] + specbuf[i]

  6. Return bytes encoded (from LSB)
```

### 8. File placement

All new code goes in `src/minidsp_steg.c` alongside the existing LSB and BFSK implementations. No new source files needed — the spectrogram text and resampler are already compiled into `libminidsp.a`.

## Risks / Trade-offs

**[Output sample rate differs from input]** The spectext method always outputs 48 kHz, even if the host was 44.1 kHz. This is a behavior change from LSB/BFSK which preserve the input rate. Documented clearly in the API and guide.

**[Amplitude headroom]** Adding 0.02-amplitude ultrasonic tones to a near-full-scale host could cause minor clipping at constructive peaks. Mitigated by the very low amplitude (worst case: 1.02 peak, negligible distortion). Could add post-mix soft clipping if needed, but likely unnecessary.

**[Memory allocation during encode]** The resampler allocates a polyphase table (~256 KB), and the spectrogram text buffer requires `out_len` doubles. For a 60-second 48 kHz signal, that's ~2.3 MB for the spec buffer plus ~2.3 MB for the resampled host. Acceptable for offline use.

**[Detection threshold tuning]** The ultrasonic energy threshold for distinguishing spectext from plain LSB needs empirical tuning. Too low: false positives on audio with ultrasonic content (rare). Too high: misses quiet spectrogram text. The 0.02 amplitude provides a clear signal well above natural noise floors.
