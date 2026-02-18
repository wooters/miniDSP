# MiniDSP

[![Docs](https://github.com/wooters/miniDSP/actions/workflows/docs.yml/badge.svg)](https://github.com/wooters/miniDSP/actions/workflows/docs.yml)
[![Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://wooters.github.io/miniDSP/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://github.com/wooters/miniDSP/blob/main/LICENSE)
[![Language: C](https://img.shields.io/badge/language-C23-orange)](https://en.wikipedia.org/wiki/C23_(C_standard_revision))

A small C library of DSP (Digital Signal Processing) routines for audio applications.

**[Read the full documentation](https://wooters.github.io/miniDSP/)** -- API reference, tutorials, and interactive examples.

## What's in the box?

### Signal Processing (minidsp.h)

- **GCC-PHAT** -- estimate the time delay between two microphone signals using Generalized Cross-Correlation with Phase Transform.  This is the core of acoustic source localisation.
- **Magnitude spectrum** -- compute |X(k)| from a real signal using the FFT; the foundation of frequency-domain analysis.
- **Power spectral density** -- compute |X(k)|^2 / N (periodogram); shows how signal power distributes across frequencies.
- **Phase spectrum** -- compute arg(X(k)) in radians; reveals the timing of each frequency component and is a prerequisite for phase-vocoder effects.
- **Spectrogram (STFT)** -- sliding-window FFT producing a time-frequency magnitude matrix; the standard tool for visualising time-varying audio.
- **Signal measurements** -- energy, power, power in dB, normalised entropy.
- **Scaling & AGC** -- linear range mapping, automatic gain control.
- **Hanning window** -- smooth windowing function for FFT analysis.
- **Sine wave generator** -- pure tone at a given frequency and amplitude; the "hello world" of DSP.
- **White noise generator** -- Gaussian random samples with configurable seed; used to test filters and measure impulse responses.
- **Impulse generator** -- single unit-sample spike at a given position; reveals a system's impulse response directly.
- **Chirp generators** -- linear and logarithmic frequency sweeps; great for testing filter magnitude response across a frequency range.
- **Square wave generator** -- bipolar square wave at a given frequency; demonstrates odd harmonics and Gibbs phenomenon.
- **Sawtooth wave generator** -- linear ramp waveform at a given frequency; contains both odd and even harmonics.

### Biquad Filters (biquad.h)

Seven classic audio filter types, all based on Robert Bristow-Johnson's [Audio EQ Cookbook](https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt):

- Low-pass, High-pass, Band-pass, Notch
- Peaking EQ, Low shelf, High shelf

### File I/O (fileio.h)

- Read audio files in any format supported by libsndfile (WAV, FLAC, AIFF, OGG, etc.)
- Write audio to WAV (IEEE float for lossless DSP round-trips)
- Write feature vectors in NumPy `.npy` format (for Python interop)
- Write feature vectors in safetensors format (for ML pipelines)
- Write feature vectors in HTK binary format (deprecated)

### Live Audio I/O (liveio.h)

- Record from the microphone and play back to speakers via PortAudio
- Non-blocking API with callback support


## Quick examples

### Detect the delay between two signals

```c
#include "minidsp.h"

/* Two 4096-sample signals captured by spatially separated microphones */
double mic_a[4096], mic_b[4096];

/* Estimate the delay in samples (+/- 50 sample search window) */
int delay = MD_get_delay(mic_a, mic_b, 4096, NULL, 50, PHAT);

printf("Signal B is %d samples behind signal A\n", delay);

/* Clean up FFTW resources when done */
MD_shutdown();
```

### Compute the magnitude spectrum

```c
#include "minidsp.h"

double signal[1024];
// ... fill signal with audio samples ...

unsigned num_bins = 1024 / 2 + 1;  /* 513 unique frequency bins */
double *mag = malloc(num_bins * sizeof(double));
MD_magnitude_spectrum(signal, 1024, mag);

/* mag[k] = |X(k)|, where frequency = k * sample_rate / 1024 */

free(mag);
MD_shutdown();
```

A full example with Hanning windowing is in `examples/magnitude_spectrum.c`.
Run it to generate an interactive HTML plot (Plotly.js + D3.js):

```sh
make -C examples plot
open examples/magnitude_spectrum.html    # interactive: zoom, pan, hover for values
```

For a step-by-step walkthrough of the DSP concepts, see the
[Magnitude Spectrum tutorial](https://wooters.github.io/miniDSP/magnitude-spectrum.html).

### Compute the power spectral density

```c
#include "minidsp.h"

double signal[1024];
// ... fill signal with audio samples ...

unsigned num_bins = 1024 / 2 + 1;  /* 513 unique frequency bins */
double *psd = malloc(num_bins * sizeof(double));
MD_power_spectral_density(signal, 1024, psd);

/* psd[k] = |X(k)|^2 / N  (power at frequency k * sample_rate / 1024) */

free(psd);
MD_shutdown();
```

A full example with Hanning windowing and one-sided PSD conversion is in
`examples/power_spectral_density.c`.
See the [PSD tutorial](https://wooters.github.io/miniDSP/power-spectral-density.html)
for a detailed explanation.

### Compute a spectrogram (STFT)

```c
#include "minidsp.h"

double signal[32000];
// ... fill signal with 2 s of audio at 16 kHz ...

unsigned N   = 512;   /* 32 ms window */
unsigned hop = 128;   /* 8 ms hop (75% overlap) */

unsigned num_frames = MD_stft_num_frames(32000, N, hop);  /* 247 */
unsigned num_bins   = N / 2 + 1;                          /* 257 */

double *mag = malloc(num_frames * num_bins * sizeof(double));
MD_stft(signal, 32000, N, hop, mag);

/* mag[f * num_bins + k] = |X_f(k)|
 * Time of frame f:   time_s  = (double)(f * hop) / 16000.0
 * Frequency of bin k: freq_hz = (double)k * 16000.0 / N       */

free(mag);
MD_shutdown();
```

A full example generating an interactive HTML heatmap is in
`examples/spectrogram.c`.
See the [Spectrogram tutorial](https://wooters.github.io/miniDSP/stft-spectrogram.html)
for a step-by-step explanation.

### Filter audio with a low-pass biquad

```c
#include "biquad.h"

/* Create a 1 kHz low-pass filter at 44.1 kHz sample rate, 1-octave bandwidth */
biquad *lpf = BiQuad_new(LPF, 0.0, 1000.0, 44100.0, 1.0);

/* Process each audio sample */
for (int i = 0; i < num_samples; i++) {
    output[i] = BiQuad(input[i], lpf);
}

free(lpf);
```

## Build and Test

### Dependencies

Install the following libraries before building:

| Library | Purpose | Debian/Ubuntu | macOS (Homebrew) |
|---------|---------|--------------|------------------|
| [FFTW3](http://www.fftw.org/) | Fast Fourier Transform | `apt install libfftw3-dev` | `brew install fftw` |
| [PortAudio](http://portaudio.com/) | Live audio I/O | `apt install portaudio19-dev` | `brew install portaudio` |
| [libsndfile](http://libsndfile.github.io/libsndfile/) | Audio file reading | `apt install libsndfile1-dev` | `brew install libsndfile` |
| [Doxygen](https://www.doxygen.nl/) | API docs generation (optional) | `apt install doxygen` | `brew install doxygen` |
| [Apple container](https://github.com/apple/container) | Linux container testing (optional) | — | macOS 26+ built-in |

The Makefiles auto-detect Homebrew paths on macOS (both Apple Silicon and Intel).

On Ubuntu, GCC 14 or later is required for `-std=c23` support. Ubuntu 24.04 ships GCC 13 by default, so install `gcc-14` explicitly (`apt install gcc-14`).

### Compile the library

```sh
make            # builds libminidsp.a
```

### Run the test suite

```sh
make test       # builds and runs all tests
```

### Test inside an Ubuntu container

To verify the library builds and passes all tests on Linux (Ubuntu 24.04 with GCC 14):

```sh
make container-test   # builds image, then runs make test inside the container
```

This requires the Apple [container](https://github.com/apple/container) CLI on macOS 26+.

### Generate API documentation

```sh
make docs       # generates HTML docs in docs/html
```

### Install git hooks

A pre-push hook is included that runs `make test` and `make container-test` before allowing pushes to `main`:

```sh
make install-hooks
```

## Roadmap

See [TODO.md](TODO.md) for planned features -- FFT spectrum analysis, signal generators, FIR filters, window functions, simple effects, pitch detection, and more.

