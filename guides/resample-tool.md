# resample -- Sample Rate Converter {#resample-tool}

**resample** converts a mono audio file from one sample rate to another using
the miniDSP polyphase sinc resampler.  It reads any format supported by
libsndfile (WAV, FLAC, AIFF, OGG, etc.) and writes WAV output (IEEE float).

## Build

From the repository root:

```sh
make tools
```

This compiles `tools/resample/resample` (requires FFTW3 and libsndfile).

## Usage

```
resample [-z N] [-b F] <input> <rate> <output.wav>
```

| Argument | Description |
|----------|-------------|
| `<input>` | Input audio file (must be mono) |
| `<rate>` | Target sample rate in Hz (positive integer) |
| `<output.wav>` | Output file path (WAV format only) |

## Options

| Flag | Default | What it does |
|------|---------|--------------|
| `-z N` | 32 | Number of sinc zero-crossings per side.  More zero-crossings produce a sharper cutoff and better stopband rejection at the cost of speed. |
| `-b F` | 10.0 | Kaiser window beta parameter.  Higher values widen the mainlobe but deepen stopband attenuation (10.0 gives >100 dB rejection). |

The defaults work well for most use cases.  Bump them up (e.g. `-z 64 -b 14.0`)
for mastering-quality conversion where you want even cleaner anti-aliasing.

## Examples

**Upsample from 16 kHz to 48 kHz:**

```sh
./tools/resample/resample speech.wav 48000 speech_48k.wav
```

**Downsample to 8 kHz (telephone quality):**

```sh
./tools/resample/resample recording.wav 8000 recording_8k.wav
```

The resampler automatically applies an anti-aliasing lowpass filter when
downsampling — no separate filtering step is needed.

**High-quality conversion with custom parameters:**

```sh
./tools/resample/resample -z 64 -b 14.0 master.wav 96000 master_96k.wav
```

## How it works

```
WAV/FLAC/... --> FIO_read_audio() --> float*
                                        |
                              float-to-double conversion
                                        |
                              MD_resample() (polyphase sinc)
                                        |
                              double-to-float conversion
                                        |
                              FIO_write_wav() --> output.wav
```

The tool reads the input file via libsndfile, converts the samples to
double precision for the DSP pipeline, resamples using `MD_resample()`,
converts back to float, and writes the output as IEEE float WAV.

Buffers are freed eagerly after each stage to keep peak memory at roughly
2x the larger of the input or output signal.

**Same-rate detection:** if the input is already at the target rate, the tool
copies the audio directly and prints a note — no resampling overhead.

**Constraints:**
- Input must be mono (single-channel).  Use `sox` to split multi-channel files.
- Output is always WAV format.  A warning is printed if the output filename
  does not end in `.wav`.

## The resampling core

The core algorithm is a 512-phase polyphase sinc interpolation filter with
Kaiser windowing.  For a detailed explanation of the math, including the
Bessel function, normalised sinc, and filter design, see the
\ref resampling "Sample Rate Conversion" tutorial.

\snippet resample.c resample-core
