# Audio Steganography {#audio-steganography}

[Steganography](https://en.wikipedia.org/wiki/Steganography) is the practice
of hiding a secret message inside an innocuous-looking cover medium.
**Audio steganography** hides data inside an audio signal so that a casual
listener hears only the original sound, while a decoder can extract the
hidden payload.

miniDSP provides three complementary methods in `src/minidsp_steg.c`,
demonstrated in `tools/audio_steg/audio_steg.c`:

| Method | Identifier | Capacity | Robustness | Audibility |
|:-------|:-----------|:---------|:-----------|:-----------|
| **LSB** (Least Significant Bit) | `MD_STEG_LSB` | High (~1 bit/sample) | Fragile | Inaudible (~-90 dB) |
| **Frequency-band** (BFSK) | `MD_STEG_FREQ_BAND` | Lower (~2.6 kbit/s) | Moderate | Near-inaudible (ultrasonic) |
| **Spectrogram text** (hybrid) | `MD_STEG_SPECTEXT` | ~4 chars/sec visual | Fragile (like LSB) | Inaudible (ultrasonic) |

Build and run the self-test from the repository root:

```sh
make -C tools/audio_steg
cd tools/audio_steg && ./audio_steg
```

---

## Message framing

Both methods prepend a **32-bit little-endian header** before the payload.
Bits 0–30 hold the message byte count; bit 31 is a **payload type flag**
(0 = text, 1 = binary).  This allows the decoder to recover the message
without knowing its length in advance, and enables `MD_steg_detect()` to
identify the payload type:

```
[ bit 31: type flag | bits 0-30: msg_len (LE) ] [ 8 * msg_len bits: payload ]
```

Each bit of the header and payload is encoded independently using the
chosen method.  Bits within each byte are transmitted LSB-first.

---

## Method 1: Least Significant Bit (LSB)

### The idea

Audio samples are typically stored as 16-bit integers (-32768 to +32767).
The least significant bit of each sample contributes only ±1 to a range of
65536 — a change of about -90 dB relative to full scale.  By replacing the
LSB of each sample with a message bit, we embed data that is completely
inaudible.

### Signal model

The host signal \f$x[n] \in [-1, 1]\f$ is quantised to 16-bit PCM:

\f[
p[n] = \mathrm{round}(x[n] \times 32767)
\f]

The LSB is then overwritten with message bit \f$b_k\f$:

\f[
p'[n] = (p[n] \mathbin{\&} \sim 1) \mathbin{|} b_k
\f]

and the stego sample is converted back to double:

\f[
y[n] = p'[n] \;/\; 32767
\f]

The maximum distortion per sample is:

\f[
|y[n] - x[n]| \leq \frac{1}{32767} \approx 3.05 \times 10^{-5}
\f]

**Reading the formula in C:**

```c
// x[n] -> host[i],  p[n] -> pcm,  b_k -> bit,  y[n] -> output[i]
int pcm = (int)(host[i] * 32767.0);  // quantise to 16-bit
pcm = (pcm & ~1) | bit;              // overwrite LSB
output[i] = (double)pcm / 32767.0;   // convert back
```

### Capacity

One bit per sample, minus the 32-bit header:

\f[
C_{\text{LSB}} = \frac{N - 32}{8} \text{ bytes}
\f]

where \f$N\f$ is the signal length in samples.

**Reading the formula in C:**

```c
// N -> signal_len,  C_LSB -> capacity
unsigned capacity = (signal_len - 32) / 8;
```

For a 3-second signal at 44.1 kHz (\f$N = 132300\f$):

\f[
C_{\text{LSB}} = \frac{132300 - 32}{8} = 16533 \text{ bytes} \approx 16 \text{ KB}
\f]

**LSB capacity at 44.1 kHz by audio duration:**

| Audio duration | Samples | Max payload | Equivalent |
|:---------------|--------:|------------:|:-----------|
| 1 s            | 44 100  | 5 508 B     | ~5 KB (small config file) |
| 5 s            | 220 500 | 27 558 B    | ~27 KB (web thumbnail) |
| 30 s           | 1 323 000 | 165 371 B | ~161 KB (high-res photo) |
| 1 min          | 2 646 000 | 330 746 B | ~323 KB (short PDF) |
| 5 min          | 13 230 000 | 1 653 746 B | ~1.6 MB (multi-page document) |
| 10 min         | 26 460 000 | 3 307 496 B | ~3.2 MB (zip archive) |
| 30 min         | 79 380 000 | 9 922 496 B | ~9.5 MB (high-res image set) |
| 1 hour         | 158 760 000 | 19 844 996 B | ~18.9 MB (small software package) |

### Listening comparison

**Original host signal** (440 Hz sine, 3 seconds):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="steg_host.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**After LSB encoding** (message hidden inside):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="steg_lsb.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

The two are perceptually identical.  The difference signal (host minus stego)
is pure quantisation noise at -90 dB:

\htmlonly
<iframe src="steg_lsb_diff.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

### Trade-offs

| Advantage | Disadvantage |
|:----------|:-------------|
| Very high capacity | Destroyed by any lossy compression (MP3, AAC, Opus) |
| Zero audible distortion | Destroyed by resampling or sample-rate conversion |
| Simple, fast implementation | Destroyed by amplitude scaling or normalisation |
| Works at any sample rate | Requires lossless transport (WAV, FLAC) |

---

## Method 2: Frequency-Band Modulation (BFSK)

### The idea

Human hearing sensitivity falls off sharply above ~16 kHz, and most adults
cannot hear tones above 18 kHz.  By adding low-amplitude tones in the
18–20 kHz "near-ultrasonic" band, we can encode data that is effectively
inaudible.

The encoding uses
[Binary Frequency-Shift Keying (BFSK)](https://en.wikipedia.org/wiki/Frequency-shift_keying):
each bit is represented by a short burst ("chip") of a sinusoidal tone at
one of two carrier frequencies.

### Carrier frequencies

| Bit value | Carrier frequency |
|:---------:|:-----------------:|
| 0         | 18500 Hz          |
| 1         | 19500 Hz          |

Both carriers are above the typical hearing threshold, and the 1 kHz
separation provides reliable discrimination during decoding.

### Chip duration

Each bit occupies a **3 ms chip** — a burst of \f$C\f$ samples:

\f[
C = \left\lfloor \frac{3.0 \times f_s}{1000} \right\rfloor
\f]

**Reading the formula in C:**

```c
// C -> chip_samples,  fs -> sample_rate
unsigned chip_samples = (unsigned)(3.0 * sample_rate / 1000.0);
```

At 44.1 kHz, \f$C = 132\f$ samples per chip.

### Encoding

For each bit \f$b_k\f$, a sine burst at the selected carrier frequency
is added to the host signal at amplitude \f$A = 0.02\f$ (-34 dB):

\f[
y[n] = x[n] + A \sin\!\bigl(2\pi\, f_{b_k}\, (n - n_0) / f_s\bigr),
\qquad n \in [n_0,\; n_0 + C)
\f]

where \f$n_0 = k \cdot C\f$ is the start sample of chip \f$k\f$ and
\f$f_{b_k}\f$ is 18500 Hz (bit 0) or 19500 Hz (bit 1).

**Reading the formula in C:**

```c
// A -> TONE_AMP (0.02),  f_bk -> freq,  n0 -> start,  fs -> sample_rate
// x[n] -> output[start+s] (already contains host), y[n] -> output[start+s]
for (unsigned s = 0; s < chip_samples; s++) {
    double t = (double)s / sample_rate;
    output[start + s] += 0.02 * sin(2.0 * M_PI * freq * t);
}
```

### Decoding

Each chip is correlated against both carrier frequencies.  The carrier
with the larger absolute correlation determines the bit value:

\f[
r_f = \sum_{s=0}^{C-1} y[n_0 + s] \,\sin\!\bigl(2\pi\, f \, s / f_s\bigr)
\f]

\f[
b_k = \begin{cases} 1 & |r_{19500}| > |r_{18500}| \\ 0 & \text{otherwise} \end{cases}
\f]

**Reading the formula in C:**

```c
// r_f -> corr_lo / corr_hi,  y[n0+s] -> stego[start+s]
double corr_lo = 0.0, corr_hi = 0.0;
for (unsigned s = 0; s < chip_samples; s++) {
    double t = (double)s / sample_rate;
    corr_lo += stego[start + s] * sin(2.0 * M_PI * 18500.0 * t);
    corr_hi += stego[start + s] * sin(2.0 * M_PI * 19500.0 * t);
}
unsigned bit = (fabs(corr_hi) > fabs(corr_lo)) ? 1 : 0;
```

### Capacity

\f[
C_{\text{freq}} = \frac{\lfloor N / C \rfloor - 32}{8} \text{ bytes}
\f]

**Reading the formula in C:**

```c
// N -> signal_len,  C -> chip_samples
unsigned total_chips = signal_len / chip_samples;
unsigned capacity = (total_chips - 32) / 8;
```

At 44.1 kHz with a 3-second signal (\f$N = 132300\f$, \f$C = 132\f$):

\f[
C_{\text{freq}} = \frac{\lfloor 132300 / 132 \rfloor - 32}{8}
                = \frac{1002 - 32}{8} = 121 \text{ bytes}
\f]

### Listening comparison

**After frequency-band encoding** (same host, message hidden via BFSK):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="steg_freq.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

The added carriers at 18.5/19.5 kHz are above most listeners' hearing range.

**Spectrogram** showing the hidden BFSK signal above the 440 Hz host tone:

\htmlonly
<iframe src="steg_freq_spectrogram.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

The faint horizontal bands near the top of the spectrogram are the BFSK
carriers.  The main 440 Hz tone dominates the audible range.

### Trade-offs

| Advantage | Disadvantage |
|:----------|:-------------|
| Survives mild additive noise | Lower capacity than LSB |
| Frequency-domain robustness | Requires sample_rate >= 40 kHz |
| Inaudible to most listeners | May be audible to young listeners with excellent high-frequency hearing |
| Amenable to spectral analysis | Vulnerable to low-pass filtering above 18 kHz |

---

## Method 3: Spectrogram Text (spectext)

### The idea

What if a hidden message were visible to the *human eye* as well as recoverable
by machine?  The **spectext** method combines LSB data encoding (for reliable
machine decode) with **spectrogram text art** in the 18--23.5 kHz ultrasonic
band (for visual verification).  Open the stego file in any spectrogram viewer
and the message is spelled out in the high frequencies — while a listener hears
nothing unusual.

The spectrogram art also acts as a **tamper indicator**: if the text is intact,
the LSB data likely is too.

### Encode pipeline

```
host.wav ──┐
(any SR)   │
           ▼
    ┌──────────────┐     ┌────────────────────┐
    │ MD_resample()│────▶│  host @ 48 kHz     │
    │ to 48 kHz    │     └────────┬───────────┘
    └──────────────┘              │
                                  ▼
                     ┌───────────────────────────┐
                     │ MD_lowpass_brickwall()     │
                     │ cutoff = original_SR / 2   │
                     │ (eliminates resampler      │
                     │  spectral images)          │
                     └────────────┬──────────────┘
                                  │
"SECRET" ──┐                      │
           ▼                      │
    ┌──────────────────────┐      │
    │MD_spectrogram_text() │      │
    │ freq: 18–23.5 kHz    │      │
    │ 30 ms / column       │      │
    │ amplitude: 0.02      │      │
    └──────────┬───────────┘      │
               │ mix (add)        │
               ▼                  ▼
           ┌──────────────────────────┐
           │  host + spectrogram art  │
           └────────────┬─────────────┘
                        │ LSB encode (last step)
                        ▼
                   stego.wav (48 kHz)
```

The spectrogram art is mixed into the host **before** LSB encoding, so the
LSB bits remain undisturbed.  Decode simply reads the LSB channel.

### Automatic upsampling and spectral cleanup

The spectrogram art uses the 18--23.5 kHz band, which requires a Nyquist
frequency of at least 23.5 kHz (sample rate >= 47 kHz).  If the input host
is below 48 kHz, the encoder automatically upsamples it using
`MD_resample()`.  After upsampling, `MD_lowpass_brickwall()` is applied at
the original Nyquist frequency to eliminate any residual spectral images
from the resampler's transition band.  This ensures the 18--23.5 kHz band
is completely clean before the spectrogram text is mixed in, so the hidden
message is clearly readable in any spectrogram viewer.  The output is
always 48 kHz.

### Fixed column width and capacity

Each character in the bitmap font is 8 columns wide (5 data + 3 spacing).
Each column occupies a fixed **30 ms** of audio, giving 240 ms per character:

\f[
C_{\text{spectext}} = \left\lfloor \frac{D}{0.24} \right\rfloor \text{ characters}
\f]

where \f$D\f$ is the host signal duration in seconds.

**Reading the formula in C:**

```c
// D -> duration_sec,  C_spectext -> vis_chars
double duration_sec = (double)signal_len / sample_rate;
unsigned vis_chars = (unsigned)(duration_sec / 0.24);
```

**Visual capacity by audio duration:**

| Audio duration | Max visible chars |
|:---------------|------------------:|
| 3 s            | 12                |
| 10 s           | 41                |
| 30 s           | 125               |
| 60 s           | 250               |

### Frequency mapping and amplitude

The 7 rows of the 5x7 bitmap font are mapped linearly across the
18--23.5 kHz band.  Row 0 (top of character) maps to 23.5 kHz;
row 6 (bottom) maps to 18 kHz.

The spectrogram text is generated at full amplitude by `MD_spectrogram_text()`
(normalised to 0.9 peak), then scaled to **0.02** (~-34 dB) before mixing.
This is loud enough to be clearly visible in a spectrogram but completely
inaudible — most adults cannot hear above 18 kHz.

### Visual truncation

If the message is longer than the visual capacity, the spectrogram art shows
only the first N characters that fit.  The full message is always recoverable
via the LSB data channel, which has much higher capacity (~5.5 KB/sec at
48 kHz).  For binary payloads, the spectrogram art shows `[BIN <N>B]` as a
label.

### Listening comparison

**After spectext encoding** (same host, message "miniDSP" hidden via spectext):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="steg_spectext.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

The ultrasonic tones at 18--23.5 kHz are far above the audible range.

**Spectrogram** showing "miniDSP" rendered as text art in the ultrasonic band:

\htmlonly
<iframe src="steg_spectext_spectrogram.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

The host audio — a [TIMIT](https://en.wikipedia.org/wiki/TIMIT) sentence
("Don't ask me to carry an oily rag like that.") — is visible at the bottom
of the spectrogram.  The text "miniDSP" is rendered in the 18--23.5 kHz band
near the top, using the 5x7 bitmap font from `MD_spectrogram_text()`.

### Trade-offs

| Advantage | Disadvantage |
|:----------|:-------------|
| Human-readable visual watermark | Lower visual capacity than LSB data capacity |
| Machine-readable round-trip via LSB | Requires 48 kHz output (auto-upsampled) |
| Visual tamper indicator | Destroyed by lossy compression (like LSB) |
| Completely inaudible | Destroyed by low-pass filtering above 18 kHz |

---

## Embedding binary data

The string-based API (`MD_steg_encode` / `MD_steg_decode`) uses null-terminated
C strings, which cannot represent binary data containing `0x00` bytes.  To hide
arbitrary binary payloads — images, compressed archives, cryptographic keys —
use the byte-oriented API:

- `MD_steg_encode_bytes()` — accepts a raw byte buffer and length
- `MD_steg_decode_bytes()` — returns raw bytes without null termination

### Visual demo

**Space invader** — a tiny 110x80 RGB PNG (332 bytes) hidden inside a short
440 Hz sine using LSB.  The recovered image is bit-identical to the original:

\htmlonly
<div style="display:flex; align-items:center; gap:2em; flex-wrap:wrap; margin:0.8em 0;">
  <div style="text-align:center;">
    <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Original</div>
    <img src="space_invader.png" alt="Space invader sprite" style="width:110px; image-rendering:pixelated;">
  </div>
  <div style="font-size:1.6em; color:#999;">&rarr;</div>
  <div style="text-align:center;">
    <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Recovered from audio</div>
    <img src="steg_recovered_invader.png" alt="Recovered space invader" style="width:110px; image-rendering:pixelated;">
  </div>
</div>
<div style="margin:0.5em 0;">
  <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Stego audio (image hidden inside):</div>
  <audio controls style="margin:0;">
    <source src="steg_lsb_invader.wav" type="audio/wav">
    <em>Your browser does not support the audio element.</em>
  </audio>
</div>
\endhtmlonly

**QR code** — a 165x165 grayscale PNG (486 bytes) encoding this repository's
URL, doubly encoded: data &rarr; QR &rarr; audio:

\htmlonly
<div style="display:flex; align-items:center; gap:2em; flex-wrap:wrap; margin:0.8em 0;">
  <div style="text-align:center;">
    <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Original</div>
    <img src="minidsp_qr.png" alt="QR code" style="width:132px; image-rendering:pixelated;">
  </div>
  <div style="font-size:1.6em; color:#999;">&rarr;</div>
  <div style="text-align:center;">
    <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Recovered from audio</div>
    <img src="steg_recovered_qr.png" alt="Recovered QR code" style="width:132px; image-rendering:pixelated;">
  </div>
</div>
<div style="margin:0.5em 0;">
  <div style="font-size:0.85em; color:#666; margin-bottom:0.3em;">Stego audio (QR hidden inside):</div>
  <audio controls style="margin:0;">
    <source src="steg_lsb_qr.wav" type="audio/wav">
    <em>Your browser does not support the audio element.</em>
  </audio>
</div>
\endhtmlonly

### Minimum samples

For LSB encoding, each data byte requires 8 samples, plus a 32-bit header:

\f[
N_{\min} = 8L + 32
\f]

where \f$L\f$ is the data length in bytes.

**Reading the formula in C:**

```c
// L -> data_len,  N_min -> min_samples
unsigned min_samples = data_len * 8 + 32;
```

For a 332-byte PNG image: \f$N_{\min} = 8 \times 332 + 32 = 2688\f$ samples —
well under a second at any common sample rate.

### Example: hiding an image in audio

**Space invader** — a tiny 110x80 RGB PNG (332 bytes) makes a good test payload.
It is small enough to embed in a fraction of a second of audio using LSB:

```sh
# Encode the image (auto-generates a host signal)
./audio_steg --encode-image lsb space_invader.png -o steg_invader.wav

# Decode to recover the image
./audio_steg --decode-image lsb steg_invader.wav -o recovered_invader.png

# Verify
cmp space_invader.png recovered_invader.png && echo "Identical"
```

**QR code** — a 165x165 1-bit grayscale PNG (486 bytes) containing a URL to
this repository.  This is a "double encoding": data encoded as a QR code,
then the QR code hidden inside audio:

```sh
./audio_steg --encode-image lsb minidsp_qr.png -o steg_qr.wav
./audio_steg --decode-image lsb steg_qr.wav -o recovered_qr.png
cmp minidsp_qr.png recovered_qr.png && echo "Identical"
```

**Quick example** — encode and decode a PNG in C:

```c
#include "minidsp.h"
#include <stdio.h>

// Read the image file
FILE *fp = fopen("space_invader.png", "rb");
fseek(fp, 0, SEEK_END);
unsigned len = (unsigned)ftell(fp);
fseek(fp, 0, SEEK_SET);
unsigned char *img = malloc(len);
fread(img, 1, len, fp);
fclose(fp);

// Encode into audio
unsigned N = len * 8 + 32 + 1024;  // payload + header + margin
double *host  = malloc(N * sizeof(double));
double *stego = malloc(N * sizeof(double));
MD_sine_wave(host, N, 0.8, 440.0, 44100.0);
MD_steg_encode_bytes(host, stego, N, 44100.0, img, len, MD_STEG_LSB);

// Decode
unsigned char *recovered = malloc(len);
MD_steg_decode_bytes(stego, N, 44100.0, recovered, len, MD_STEG_LSB);
// recovered[0..len-1] == img[0..len-1]
```

---

## API

### Capacity

```c
unsigned MD_steg_capacity(unsigned signal_len, double sample_rate, int method);
```

Returns the maximum number of message bytes that can be hidden.

### Encode

```c
unsigned MD_steg_encode(const double *host, double *output,
                        unsigned signal_len, double sample_rate,
                        const char *message, int method);
```

| Parameter     | Description |
|:--------------|:------------|
| `host`        | Input host signal (not modified). |
| `output`      | Output stego signal (caller-allocated, same length). |
| `signal_len`  | Number of samples. |
| `sample_rate` | Sample rate in Hz. |
| `message`     | Null-terminated secret message. |
| `method`      | `MD_STEG_LSB`, `MD_STEG_FREQ_BAND`, or `MD_STEG_SPECTEXT`. |

Returns the number of message bytes encoded (0 on failure).

### Decode

```c
unsigned MD_steg_decode(const double *stego, unsigned signal_len,
                        double sample_rate,
                        char *message_out, unsigned max_msg_len,
                        int method);
```

| Parameter     | Description |
|:--------------|:------------|
| `stego`       | The stego signal containing the hidden message. |
| `signal_len`  | Number of samples. |
| `sample_rate` | Sample rate in Hz. |
| `message_out` | Output buffer (caller-allocated, null-terminated on return). |
| `max_msg_len` | Size of buffer including null terminator. |
| `method`      | `MD_STEG_LSB`, `MD_STEG_FREQ_BAND`, or `MD_STEG_SPECTEXT`. |

Returns the number of message bytes decoded (0 if none found).

### Encode bytes

```c
unsigned MD_steg_encode_bytes(const double *host, double *output,
                              unsigned signal_len, double sample_rate,
                              const unsigned char *data, unsigned data_len,
                              int method);
```

| Parameter     | Description |
|:--------------|:------------|
| `host`        | Input host signal (not modified). |
| `output`      | Output stego signal (caller-allocated, same length). |
| `signal_len`  | Number of samples. |
| `sample_rate` | Sample rate in Hz. |
| `data`        | Pointer to the binary data to hide. |
| `data_len`    | Length of data in bytes. |
| `method`      | `MD_STEG_LSB`, `MD_STEG_FREQ_BAND`, or `MD_STEG_SPECTEXT`. |

Returns the number of data bytes encoded (0 on failure).

### Decode bytes

```c
unsigned MD_steg_decode_bytes(const double *stego, unsigned signal_len,
                              double sample_rate,
                              unsigned char *data_out, unsigned max_len,
                              int method);
```

| Parameter     | Description |
|:--------------|:------------|
| `stego`       | The stego signal containing the hidden data. |
| `signal_len`  | Number of samples. |
| `sample_rate` | Sample rate in Hz. |
| `data_out`    | Output buffer for the decoded bytes (caller-allocated). |
| `max_len`     | Maximum number of bytes to write to buffer. |
| `method`      | `MD_STEG_LSB`, `MD_STEG_FREQ_BAND`, or `MD_STEG_SPECTEXT`. |

Returns the number of data bytes decoded (0 if none found).

### Detect

```c
int MD_steg_detect(const double *signal, unsigned signal_len,
                   double sample_rate, int *payload_type_out);
```

Inspects a signal and determines which steganography method (if any) was used.
Returns `MD_STEG_LSB`, `MD_STEG_FREQ_BAND`, `MD_STEG_SPECTEXT`, or `-1` if
no hidden payload is found.  The optional `payload_type_out` receives `MD_STEG_TYPE_TEXT` (0) or
`MD_STEG_TYPE_BINARY` (1).

| Parameter          | Description |
|:-------------------|:------------|
| `signal`           | The signal to inspect. |
| `signal_len`       | Length of the signal in samples. |
| `sample_rate`      | Sample rate in Hz. |
| `payload_type_out` | If non-null, receives the payload type flag. |

**How it works:** The function probes the first 32 samples (LSB) or 32 BFSK
chips (frequency-band) to extract the header.  A header is considered valid
when the decoded length is positive and fits the signal capacity.  For BFSK,
the average correlation must also exceed a minimum threshold to avoid false
positives.  If both methods claim a valid header, BFSK wins (harder to
trigger by accident).

**Quick example:**

```c
int payload_type;
int method = MD_steg_detect(signal, signal_len, 44100.0, &payload_type);
if (method == MD_STEG_LSB)
    printf("LSB-encoded %s payload detected\n",
           payload_type == MD_STEG_TYPE_BINARY ? "binary" : "text");
```

---

## Quick example

**Encode and decode with LSB:**

```c
#include "minidsp.h"

double host[44100], stego[44100];
MD_sine_wave(host, 44100, 0.8, 440.0, 44100.0);

// Encode
unsigned n = MD_steg_encode(host, stego, 44100, 44100.0,
                            "secret message", MD_STEG_LSB);

// Decode
char recovered[256];
MD_steg_decode(stego, 44100, 44100.0, recovered, 256, MD_STEG_LSB);
printf("Hidden: %s\n", recovered);  // "secret message"
```

**Encode and decode with frequency band:**

```c
double host[132300], stego[132300];   // 3 s at 44.1 kHz
MD_sine_wave(host, 132300, 0.8, 440.0, 44100.0);

unsigned n = MD_steg_encode(host, stego, 132300, 44100.0,
                            "hidden!", MD_STEG_FREQ_BAND);

char recovered[256];
MD_steg_decode(stego, 132300, 44100.0, recovered, 256, MD_STEG_FREQ_BAND);
printf("Hidden: %s\n", recovered);  // "hidden!"
```

**Encode and decode with spectrogram text (spectext):**

```c
double host[132300];  // 3 s at 44.1 kHz
MD_sine_wave(host, 132300, 0.8, 440.0, 44100.0);

// Output at 48 kHz — compute required buffer size
unsigned out_len = MD_resample_output_len(132300, 44100.0, 48000.0);
double *stego = malloc(out_len * sizeof(double));

MD_steg_encode(host, stego, 132300, 44100.0, "miniDSP", MD_STEG_SPECTEXT);

// Decode from the 48 kHz output
char recovered[256];
MD_steg_decode(stego, out_len, 48000.0, recovered, 256, MD_STEG_SPECTEXT);
printf("Hidden: %s\n", recovered);  // "miniDSP"
// View stego in a spectrogram to see "miniDSP" in the 18-23.5 kHz band
free(stego);
```

---

## Example program

The tool `tools/audio_steg/audio_steg.c` provides a command-line program for
encoding and decoding steganographic messages and binary data in WAV files.

**Self-test** (no arguments):

\snippet audio_steg.c self-test

**Encode a message into a WAV file:**

\snippet audio_steg.c encode-wav

**Decode a message from a WAV file:**

\snippet audio_steg.c decode-wav

**Encode a binary file (e.g. image) into a WAV file:**

\snippet audio_steg.c encode-image

**Decode a binary file from a WAV file:**

\snippet audio_steg.c decode-image

**Auto-detect and decode** (no method needed):

\snippet audio_steg.c auto-detect

**Usage:**

```sh
# Self-test (encode + decode with both methods, verify round-trip)
./audio_steg

# Encode a text message using LSB into a default 440 Hz host
./audio_steg --encode lsb "my secret" -o stego.wav

# Encode using frequency-band into an existing WAV host
./audio_steg --encode freq "hidden" -i music.wav -o stego.wav

# Auto-detect and decode (just pass the file)
./audio_steg stego.wav

# Decode with explicit method
./audio_steg --decode lsb stego.wav
./audio_steg --decode freq stego.wav

# Decode without specifying method (auto-detect)
./audio_steg --decode stego.wav

# Encode using spectext (hybrid LSB + spectrogram art)
./audio_steg --encode spectext "miniDSP" -i music.wav -o stego.wav

# Encode a binary file (image) using LSB
./audio_steg --encode-image lsb space_invader.png -o steg_invader.wav

# Decode a binary file (auto-detect method)
./audio_steg --decode-image steg_invader.wav -o recovered.png
```

---

## Choosing a method

| Criterion | LSB | Frequency-band | Spectext |
|:----------|:----|:---------------|:---------|
| **Message size** | Up to ~16 KB/s of audio | Up to ~121 B per 3 s | ~4 chars/s (visual); LSB capacity for data |
| **Audio quality** | Imperceptible (-90 dB) | Near-imperceptible (-34 dB) | Imperceptible (ultrasonic, -34 dB) |
| **Survives lossy compression** | No | No (but tolerates noise) | No |
| **Survives additive noise** | No (bit errors) | Yes (mild noise) | No (LSB channel) |
| **Sample rate requirement** | Any | >= 40 kHz | Output always 48 kHz |
| **Visual verification** | No | No (spectrogram shows carriers) | Yes — text readable in spectrogram |
| **Best for** | Lossless pipelines (WAV/FLAC) | Light interference environments | Visual watermarking + machine decode |

For maximum capacity and fidelity in lossless pipelines, use **LSB**.
For slightly more robust hiding in near-ultrasonic bands, use
**frequency-band**.  For a human-readable visual watermark with machine-readable
data recovery, use **spectext**.
