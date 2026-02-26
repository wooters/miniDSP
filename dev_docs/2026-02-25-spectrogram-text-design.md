# Spectrogram Text Art — Design Overview

## Idea

Generate an audio signal that displays readable text when viewed as a spectrogram. Given a short string like "HELLO", synthesize sine waves whose frequencies and timing trace out the characters. The result is a WAV file that sounds like a buzzy chord but reveals the hidden message visually.

## Motivation

This feature demonstrates the duality between time-domain audio and frequency-domain visualization. It combines signal generation (already a strength of miniDSP) with the existing STFT spectrogram infrastructure to produce a striking visual result — useful for teaching, demos, and fun.

## API

One new public function in `src/minidsp_spectext.c`:

```c
unsigned MD_spectrogram_text(double *output, unsigned max_len,
                             const char *text,
                             double freq_lo, double freq_hi,
                             double duration_sec, double sample_rate);
```

The caller provides a buffer and parameters; the function fills the buffer with synthesized audio and returns the number of samples written. A built-in 5x7 bitmap font rasterizes printable ASCII characters. Each bitmap column becomes a time slice; each "on" pixel becomes a sine wave at the corresponding frequency. A 3 ms raised-cosine crossfade suppresses clicking at column boundaries. The output is normalized to 0.9 peak amplitude.

## Example Program

`examples/spectrogram_text.c` generates a WAV file, an HTML Plotly spectrogram, and a PNG image. The PNG uses either a Viridis or grayscale colormap, selected via `--colormap`. The STFT uses a large FFT and small hop size for a sharp spectrogram image (computation time is not a concern). The single-header library `stb_image_write.h` handles PNG encoding with no additional link-time dependencies.

## Implementation Plan

See `2026-02-25-spectrogram-text-plan.md`.
