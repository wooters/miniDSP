# MiniDSP

A small C library of DSP (Digital Signal Processing) routines for audio applications.

## What's in the box?

### Signal Processing (`minidsp.h`)
- **GCC-PHAT** -- estimate the time delay between two microphone signals using Generalized Cross-Correlation with Phase Transform.  This is the core of acoustic source localisation.
- **Signal measurements** -- energy, power, power in dB, normalised entropy.
- **Scaling & AGC** -- linear range mapping, automatic gain control.
- **Hanning window** -- smooth windowing function for FFT analysis.

### Biquad Filters (`biquad.h`)
Seven classic audio filter types, all based on Robert Bristow-Johnson's Audio EQ Cookbook:
- Low-pass, High-pass, Band-pass, Notch
- Peaking EQ, Low shelf, High shelf

### File I/O (`fileio.h`)
- Read audio files in any format supported by libsndfile (WAV, FLAC, AIFF, OGG, etc.)
- Write feature vectors in HTK binary format (for speech recognition pipelines)

### Live Audio I/O (`liveio.h`)
- Record from the microphone and play back to speakers via PortAudio
- Non-blocking API with callback support

## Building

### Dependencies

| Library | Purpose | Debian/Ubuntu | macOS (Homebrew) |
|---------|---------|--------------|------------------|
| [FFTW3](http://www.fftw.org/) | Fast Fourier Transform | `apt install libfftw3-dev` | `brew install fftw` |
| [PortAudio](http://portaudio.com/) | Live audio I/O | `apt install portaudio19-dev` | `brew install portaudio` |
| [libsndfile](http://libsndfile.github.io/libsndfile/) | Audio file reading | `apt install libsndfile1-dev` | `brew install libsndfile` |

The Makefiles auto-detect Homebrew paths on macOS (both Apple Silicon and Intel).

### Compile the library

```sh
make            # builds libminidsp.a
```

### Run the test suite

```sh
make test       # builds and runs all 47 tests
```

## Quick example: detect the delay between two signals

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

## Quick example: filter audio with a low-pass biquad

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

## Test suite

The test suite (`tests/test_minidsp.c`) covers every public function:

- **Dot product** -- orthogonal vectors, known values, self-dot
- **Energy / Power / dB** -- known signals, sine wave power, dB floor
- **Scaling** -- endpoints, midpoint, vector scaling, fit-within-range
- **AGC** -- target dB level achievement
- **Entropy** -- uniform, spike, zero, clip/no-clip modes
- **Hanning window** -- endpoints, peak, symmetry, range
- **GCC-PHAT** -- positive/negative/zero delays, SIMP vs PHAT weighting, multi-signal delays, FFT plan caching
- **Biquad filters** -- LPF, HPF, BPF, Notch, PEQ, Low shelf, High shelf, DC rejection

## License

MIT License.  See [LICENSE](LICENSE) for details.

## Author

Chuck Wooters -- <wooters@hey.com>
