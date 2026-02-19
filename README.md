# MiniDSP

[![Docs](https://github.com/wooters/miniDSP/actions/workflows/docs.yml/badge.svg)](https://github.com/wooters/miniDSP/actions/workflows/docs.yml)
[![Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://wooters.github.io/miniDSP/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](https://github.com/wooters/miniDSP/blob/main/LICENSE)
[![Language: C](https://img.shields.io/badge/language-C23-orange)](https://en.wikipedia.org/wiki/C23_(C_standard_revision))

A small C library of DSP (Digital Signal Processing) routines for audio applications.

**[Read the full documentation](https://wooters.github.io/miniDSP/)** -- API reference, tutorials, and interactive examples.

## What's in the box?

### Spectral Analysis (minidsp.h)

- **Magnitude spectrum** -- compute |X(k)| from a real signal using the FFT; the foundation of frequency-domain analysis.
- **Power spectral density** -- compute |X(k)|^2 / N (periodogram); shows how signal power distributes across frequencies.
- **Phase spectrum** -- compute arg(X(k)) in radians; reveals the timing of each frequency component and is a prerequisite for phase-vocoder effects.
- **Spectrogram (STFT)** -- sliding-window FFT producing a time-frequency magnitude matrix; the standard tool for visualising time-varying audio.
- **Mel filterbank** -- triangular filters spaced on the mel scale for perceptually motivated spectral features.
- **MFCCs** -- mel-frequency cepstral coefficients from log mel energies via DCT-II (C0 included).
- **Window functions** -- Hanning, Hamming, Blackman, and rectangular windows for FFT analysis trade-off studies.

### Signal Generators (minidsp.h)

- **Sine wave generator** -- pure tone at a given frequency and amplitude; the "hello world" of DSP.
- **White noise generator** -- Gaussian random samples with configurable seed; used to test filters and measure impulse responses.
- **Impulse generator** -- single unit-sample spike at a given position; reveals a system's impulse response directly.
- **Chirp generators** -- linear and logarithmic frequency sweeps; great for testing filter magnitude response across a frequency range.
- **Square wave generator** -- bipolar square wave at a given frequency; demonstrates odd harmonics and Gibbs phenomenon.
- **Sawtooth wave generator** -- linear ramp waveform at a given frequency; contains both odd and even harmonics.

### Filters (biquad.h)

Seven classic audio filter types, all based on Robert Bristow-Johnson's [Audio EQ Cookbook](https://webaudio.github.io/Audio-EQ-Cookbook/Audio-EQ-Cookbook.txt):

- Low-pass, High-pass, Band-pass, Notch
- Peaking EQ, Low shelf, High shelf

### Delay Estimation (minidsp.h)

- **GCC-PHAT** -- estimate the time delay between two microphone signals using Generalized Cross-Correlation with Phase Transform.  This is the core of acoustic source localisation.

### Signal Analysis (minidsp.h)

- **RMS** -- root mean square, the standard signal loudness measure.
- **Zero-crossing rate** -- fraction of adjacent samples that change sign; simple proxy for pitch and noisiness.
- **Autocorrelation** -- normalised self-similarity at different lags; foundation of pitch detection.
- **Peak detection** -- find local maxima above a threshold with minimum-distance suppression.
- **F0 estimation (autocorrelation)** -- estimate pitch by finding the dominant autocorrelation lag in a frequency range.
- **F0 estimation (FFT peak-pick)** -- estimate pitch from the dominant Hann-windowed spectral peak.
- **Signal mixing** -- weighted sum of two signals; needed for any multi-source demo.

### Simple Effects (minidsp.h)

- **Delay line / echo** -- circular-buffer delay with feedback; the building block of many audio effects.
- **Tremolo** -- amplitude modulation by a low-frequency oscillator.
- **Comb-filter reverb** -- feedback comb filter introducing a reverb-like decay tail.

### FIR Filters / Convolution (minidsp.h)

- **Time-domain convolution** -- direct full linear convolution for teaching and validation.
- **Moving-average filter** -- simplest causal FIR low-pass with zero-padded startup.
- **General FIR filter** -- apply arbitrary tap coefficients to build custom FIR responses.
- **FFT overlap-add convolution** -- fast full convolution for longer kernels.

### Signal Measurement (minidsp.h)

- **Signal measurements** -- energy, power, power in dB, normalised entropy.
- **Scaling & AGC** -- linear range mapping, automatic gain control.

### File I/O (fileio.h)

- Read audio files in any format supported by libsndfile (WAV, FLAC, AIFF, OGG, etc.)
- Write audio to WAV (IEEE float for lossless DSP round-trips)
- Write feature vectors in NumPy `.npy` format (for Python interop)
- Write feature vectors in safetensors format (for ML pipelines)
- Write feature vectors in HTK binary format (deprecated)

### Live Audio I/O (liveio.h)

- Record from the microphone and play back to speakers via PortAudio
- Non-blocking API with callback support

## Use in your project

Install the [dependencies](#dependencies) listed below, then clone and build:

```sh
git clone https://github.com/wooters/miniDSP.git
cd miniDSP
make
```

This produces `libminidsp.a` in the repo root.

**A minimal program** -- generate a sine wave and find its peak frequency bin:

```c
#include "minidsp.h"
#include <stdio.h>

int main(void) {
    double signal[1024];
    MD_sine_wave(signal, 1024, 1.0, 440.0, 16000.0);

    double mag[1024 / 2 + 1];
    MD_magnitude_spectrum(signal, 1024, mag);

    unsigned peak = 0;
    for (unsigned k = 1; k < 1024 / 2 + 1; k++)
        if (mag[k] > mag[peak]) peak = k;
    printf("Peak bin: %u (%.1f Hz)\n", peak, peak * 16000.0 / 1024);

    MD_shutdown();
}
```

Save the code above as `my_program.c`.

**Compile it directly:**

```sh
gcc -std=c23 -Ipath/to/miniDSP/include my_program.c \
    -Lpath/to/miniDSP -lminidsp -lfftw3 -lm -o my_program
```

**Or use a Makefile** (adapts to Homebrew on macOS automatically):

```makefile
MINIDSP_DIR = path/to/miniDSP

CC      = gcc
CFLAGS  = -std=c23 -Wall -Wextra -I$(MINIDSP_DIR)/include
LDFLAGS = -L$(MINIDSP_DIR)
LDLIBS  = -lminidsp -lfftw3 -lm

BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
  CFLAGS  += -I$(BREW_PREFIX)/include
  LDFLAGS += -L$(BREW_PREFIX)/lib
endif

my_program: my_program.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< $(LDLIBS) -o $@
```

If you use `fileio.h` for reading or writing audio files, add `-lsndfile` to `LDLIBS`.

## Quick examples

For step-by-step walkthroughs of these and other topics, see the
[Tutorials](https://wooters.github.io/miniDSP/tutorials.html) in the full documentation.

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

### Estimate fundamental frequency (F0)

```c
#include "minidsp.h"

double frame[1024];
// ... fill frame with one analysis frame of audio ...

double f0_acf = MD_f0_autocorrelation(frame, 1024, 16000.0, 80.0, 400.0);
double f0_fft = MD_f0_fft(frame, 1024, 16000.0, 80.0, 400.0);

/* Either value may be 0.0 if no reliable F0 is found. */
MD_shutdown();
```

A runnable frame-tracking example is in `examples/pitch_detection.c`.
See the [Pitch Detection tutorial](https://wooters.github.io/miniDSP/pitch-detection.html)
for method comparison and visuals.

### Compute mel energies and MFCCs

```c
#include "minidsp.h"

double frame[1024];
// ... fill frame with one analysis frame ...

double mel[26];
double mfcc[13];

MD_mel_energies(frame, 1024, 16000.0, 26, 80.0, 7600.0, mel);
MD_mfcc(frame, 1024, 16000.0, 26, 13, 80.0, 7600.0, mfcc);

/* mfcc[0] is C0 */
MD_shutdown();
```

A runnable example is in `examples/mel_mfcc.c`.
See the [Mel/MFCC tutorial](https://wooters.github.io/miniDSP/mel-mfcc.html).

### FIR filtering and convolution

```c
#include "minidsp.h"

double x[1024];      // input signal
double h[] = {0.2, 0.6, 0.2};   // simple smoothing kernel

unsigned ylen = MD_convolution_num_samples(1024, 3);
double *y_time = malloc(ylen * sizeof(double));
double *y_fft  = malloc(ylen * sizeof(double));

MD_convolution_time(x, 1024, h, 3, y_time);
MD_convolution_fft_ola(x, 1024, h, 3, y_fft);   // same output, faster for long kernels

double y_fir[1024];
MD_fir_filter(x, 1024, h, 3, y_fir);

double y_ma[1024];
MD_moving_average(x, 1024, 8, y_ma);

free(y_fft);
free(y_time);
MD_shutdown();
```

A runnable example is in `examples/fir_convolution.c`.
See the [FIR/Convolution tutorial](https://wooters.github.io/miniDSP/fir-convolution.html).

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
