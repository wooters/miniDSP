# Spectrogram Text Art {#spectrogram-text}

Generate audio that displays readable text when viewed as a spectrogram.

## The idea

A spectrogram is a time-frequency picture of a signal: the x-axis is time,
the y-axis is frequency, and brightness encodes magnitude.  If we place
sine waves at specific frequencies during specific time intervals, we can
"draw" in the spectrogram.  Given a short string like `"HELLO"`, the
library rasterises it with a tiny bitmap font, maps each pixel to a
frequency band, and synthesises the corresponding tones.  The result is a
WAV file that sounds like a buzzy chord but reveals the hidden message
visually.

## How it works

**Bitmap font** — Each printable ASCII character (32–126) is stored as a
5-column, 7-row bitmap.  One byte per column, one bit per row.  Characters
are spaced three columns apart, giving a grid width of
\f$ W = 8 \times \mathrm{len} - 3 \f$ columns and 7 rows.

**Frequency mapping** — Row 0 maps to `freq_hi` (top of the spectrogram),
row 6 maps to `freq_lo` (bottom).  Intermediate rows are linearly
interpolated:

\f[
f_r = f_\mathrm{hi} - \frac{r}{6} \, (f_\mathrm{hi} - f_\mathrm{lo})
\f]

**Reading the formula in C:**

```c
// f_hi -> freq_hi,  f_lo -> freq_lo,  r -> r,  f_r -> row_freq[r]
// The font has 7 rows, so the denominator is 6 (= 7 − 1).
double row_freq[7];
for (unsigned r = 0; r < 7; r++)
    row_freq[r] = freq_hi - (double)r / 6.0 * (freq_hi - freq_lo);
```

**Synthesis** — Each bitmap column occupies `col_samples = duration / W * sample_rate`
audio samples.  For every "on" pixel in that column a sine wave at
frequency \f$ f_r \f$ is added.  Phase is carried continuously across
columns so that sustained tones (adjacent "on" pixels in the same row)
produce a smooth sinusoid.

**Crossfade** — A 3 ms raised-cosine ramp is applied at tone onsets and
offsets (transitions between adjacent "on" and "off" columns in the same
row) to suppress broadband clicks.

**Normalisation** — After all columns are synthesised the entire signal is
scaled so that the peak absolute amplitude is 0.9.

## API

```c
unsigned MD_spectrogram_text(double *output, unsigned max_len,
                             const char *text,
                             double freq_lo, double freq_hi,
                             double duration_sec, double sample_rate);
```

**Parameters:**

| Parameter      | Description |
|----------------|-------------|
| `output`       | Caller-allocated buffer for the synthesised audio. |
| `max_len`      | Size of `output` in samples (must be >= returned value). |
| `text`         | Printable ASCII string to render (non-empty). |
| `freq_lo`      | Lowest frequency in Hz (bottom of text in spectrogram). |
| `freq_hi`      | Highest frequency in Hz (top of text in spectrogram). |
| `duration_sec` | Total signal duration in seconds. |
| `sample_rate`  | Sample rate in Hz (`freq_hi` must be < `sample_rate / 2`). |

**Returns** the number of samples written.

## Quick example

```c
#include "minidsp.h"

double buf[64000];
unsigned n = MD_spectrogram_text(buf, 64000, "HELLO",
                                 200.0, 7500.0, 2.0, 16000.0);
/* buf[0..n-1] contains the audio.  Compute its STFT to see "HELLO". */
```

## Example program

The example `examples/spectrogram_text.c` generates a WAV file, an interactive
HTML spectrogram, and a static PNG image.

**Usage:**

```sh
./spectrogram_text [TEXT] [--colormap hot|grayscale]
```

Default text is `"HELLO"`, default colormap is `hot`.

**Synthesis and STFT parameters** used by the example (so you can replicate the
spectrogram with your own tools):

| Parameter      | Value             |
|----------------|-------------------|
| Sample rate    | 16000 Hz          |
| Frequency range| 400–7300 Hz       |
| Duration       | 2.25 s            |
| Silence padding| 0.5 s before/after (default) |
| STFT FFT size  | 1024              |
| STFT hop       | 16 samples (1 ms) |
| dB range       | −80 to 0 dB       |
| Colorscale     | Hot (or grayscale) |

Outputs: `spectrogram_text.wav`, `spectrogram_text.html`, `spectrogram_text.png`.

## Viewing the result

**Listen** — the synthesised "HELLO" (buzzy chord):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="spectrogram_text_hello.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Spectrogram** — each letter is visible as a cluster of horizontal bands in
the time-frequency plane:

\htmlonly
<iframe src="spectext_hello_spectrogram.html" style="width:100%;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

**Build and run:**

```sh
make -C examples spectrogram_text
cd examples && ./spectrogram_text "HELLO"
open spectrogram_text.html     # interactive spectrogram
open spectrogram_text.png      # static image
```
