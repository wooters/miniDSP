# miniDSP

A small C library for audio DSP: signal measurement, biquad filtering, GCC-PHAT delay estimation, and audio I/O.

## Build commands

```bash
make              # Build libminidsp.a
make test         # Build and run the test suite (tests/test_minidsp)
make docs         # Generate HTML documentation (requires Doxygen)
make clean        # Remove all build artifacts (including generated docs)
```

To build and run just the test binary directly:

```bash
make -C tests test_minidsp && ./tests/test_minidsp
```

Legacy demo programs in `tests/` (gcc_phat_test, testliveio, etc.) require gnuplot_i and/or PortAudio at link time. Build them individually, e.g. `make -C tests gcc_phat_test`.

## Dependencies

- **FFTW3** (`-lfftw3`) — used by minidsp for FFT-based cross-correlation
- **libsndfile** (`<sndfile.h>`) — used by fileio for reading WAV/FLAC/AIFF
- **PortAudio** (`<portaudio.h>`, `-lportaudio`) — used by liveio for live recording/playback

The primary test suite (`test_minidsp`) links FFTW3, math, and libsndfile. PortAudio is only needed if you build the liveio-dependent programs.

## Architecture

Four modules, each a `.c`/`.h` pair:

| Module | Prefix | Purpose |
|--------|--------|---------|
| `minidsp` | `MD_` | Signal math (energy, entropy, scaling, Hanning window) and GCC-PHAT delay estimation |
| `biquad` | `BiQuad` / `BiQuad_` | Second-order IIR filter (LPF, HPF, BPF, notch, PEQ, shelving) |
| `fileio` | `FIO_` | Read/write audio via libsndfile; write feature vectors (.npy, safetensors, HTK) |
| `liveio` | `LA_` | Record/play audio via PortAudio callbacks |

Key patterns:

- **FFT plan caching** — `minidsp.c` keeps static FFTW plans and buffers. They are allocated on first use and reallocated only when signal length changes. Call `MD_shutdown()` to free them.
- **PortAudio callbacks** — `liveio` uses non-blocking PortAudio streams. `LA_record()` / `LA_play()` return immediately; audio I/O happens in background callbacks. Also supports user-supplied callbacks via `LA_record_callback()` / `LA_play_callback()`.
- **Heap-allocated filters** — `BiQuad_new()` returns a `malloc`'d struct. Caller must `free()` it.

## Code conventions

- **C23** (`-std=c23 -Wall -Wextra -pedantic`)
- All modules compile into a single static library `libminidsp.a`
- Public functions use module prefixes: `MD_`, `BiQuad_`, `FIO_`, `LA_`
- `smp_type` is `typedef double` in biquad.h (could be changed to `float`)
- Tests live in `tests/` with their own Makefile; `make test` from root delegates there
