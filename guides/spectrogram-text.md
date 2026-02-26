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
are spaced one column apart, giving a grid width of
\f$ W = 6 \times \mathrm{len} - 1 \f$ columns and 7 rows.

**Frequency mapping** — Row 0 maps to `freq_hi` (top of the spectrogram),
row 6 maps to `freq_lo` (bottom).  Intermediate rows are linearly
interpolated:

\f[
f_r = f_\mathrm{hi} - \frac{r}{6} \, (f_\mathrm{hi} - f_\mathrm{lo})
\f]

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
| `sample_rate`  | Sample rate in Hz (`freq_hi` must be <= `sample_rate / 2`). |

**Returns** the number of samples written.

## Quick example

```c
#include "minidsp.h"

double buf[64000];
unsigned n = MD_spectrogram_text(buf, 64000, "HELLO",
                                 200.0, 7500.0, 2.0, 16000.0);
/* buf[0..n-1] contains the audio.  Compute its STFT to see "HELLO". */
```

## Viewing the result

The example program `examples/spectrogram_text.c` generates three outputs:

1. **WAV file** — play back the audio to hear the buzzy chord.
2. **HTML spectrogram** — interactive Plotly heatmap; hover for exact time/frequency/dB.
3. **PNG image** — static spectrogram with Viridis or grayscale colourmap.

Build and run:

```sh
make -C examples spectrogram_text
cd examples && ./spectrogram_text "HELLO"
open spectrogram_text.html     # interactive spectrogram
open spectrogram_text.png      # static image
```
