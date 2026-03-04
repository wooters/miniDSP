# Plan: Spectrogram-to-Audio Reconstruction

## Context

Add the ability to reconstruct audio from a magnitude-only spectrogram. This requires two new library functions (inverse STFT and Griffin-Lim phase reconstruction) plus an example app that reads a miniDSP-generated PNG spectrogram and outputs a WAV file. Only our own spectrograms need to be supported (known Viridis colormap, -80 to 0 dB range, known STFT params).

## New Library Functions

### 1. `MD_istft_signal_len()` — compute output buffer size
```c
unsigned MD_istft_signal_len(unsigned num_frames, unsigned N, unsigned hop);
// Returns (num_frames - 1) * hop + N, or 0 if num_frames == 0
```

### 2. `MD_istft()` — inverse STFT with overlap-add
```c
void MD_istft(const double *mag, const double *phase,
              unsigned num_frames, unsigned N, unsigned hop,
              double *signal_out, unsigned signal_len);
```
- Takes separate magnitude + phase arrays (same layout as `MD_stft()` output: `[frame * (N/2+1) + k]`)
- FFTW c2r inverse, divide by N, Hanning synthesis window, overlap-add
- Per-sample normalization by sum of squared windows

### 3. `MD_griffin_lim()` — phase reconstruction from magnitude
```c
void MD_griffin_lim(const double *mag, unsigned num_frames,
                   unsigned N, unsigned hop, unsigned iterations,
                   double *signal_out, unsigned signal_len);
```
- Starts with zero phase, iterates: iSTFT -> forward STFT -> extract new phase -> keep original magnitude
- 32-60 iterations typical; more = better but diminishing returns

## FFTW c2r Plan Caching

New static cache in `src/minidsp_spectrum.c` (separate from existing r2c cache):
- `_istft_N`, `_istft_in` (fftw_complex), `_istft_out` (double), `_istft_plan`
- `_istft_setup(N)` / `_istft_teardown()` following existing pattern
- `_istft_teardown()` does NOT reset `_istft_N` (same rule as `_spec_teardown`)
- `MD_shutdown()` calls `_istft_teardown()` then resets `_istft_N = 0`

## Implementation Details

### Inverse STFT Algorithm
For each frame `f`:
1. Construct complex spectrum: `S(k) = mag(k) * exp(j * phase(k))`
2. FFTW c2r inverse FFT -> time-domain frame
3. Divide by N (FFTW unnormalized convention)
4. Multiply by Hanning synthesis window
5. Overlap-add into output buffer at offset `f * hop`
6. Accumulate `w[n]^2` into normalization buffer

Final step: divide each output sample by its accumulated window normalization (threshold > 1e-10 to avoid division by near-zero).

### Griffin-Lim Algorithm
1. Initialize phase to zero
2. For each iteration:
   - `MD_istft(original_mag, current_phase)` -> time signal
   - Forward STFT of time signal (using cached r2c plan + Hanning window)
   - Extract new phase via `carg()`, discard new magnitude
   - Replace current_phase with new phase
3. Final `MD_istft()` with converged phase

### Griffin-Lim Forward Analysis (inline in `minidsp_spectrum.c`)
Griffin-Lim needs forward STFT that returns **both** magnitude and phase. Since `MD_stft()` only returns magnitude, the re-analysis step is done inline using the cached `_spec_plan`, `_spec_in`, `_spec_out`, and `_stft_win` statics (all accessible within the same translation unit).

### FFTW c2r Note
`FFTW_DESTROY_INPUT` is used for the c2r plan. This is fine because `_istft_in` is repopulated fresh each frame. `FFTW_PRESERVE_INPUT` is not reliably supported for c2r transforms.

## Files to Modify

| File | Changes |
|------|---------|
| `include/minidsp.h` | Declarations + Doxygen docs for 3 new functions; update brief feature list |
| `src/minidsp_spectrum.c` | Implementations + c2r plan cache + update `MD_shutdown()` |
| `tests/test_minidsp.c` | ~10 tests (round-trip, silence, convergence, etc.) |
| `examples/Makefile` | Add `spectrogram_reconstruct` to `SNDFILE_EXAMPLES` and `plot` target |
| `.gitignore` | Add `examples/spectrogram_reconstruct{,.csv,.html,.wav,.png}` |
| `.dockerignore` | Add `examples/spectrogram_reconstruct` |
| `guides/tutorials.md` | Add `\subpage spectrogram-reconstruct` |
| `README.md` | Update feature list |

## New Files

| File | Purpose |
|------|---------|
| `third_party/stb_image.h` | Single-header PNG reader (download from nothings/stb) |
| `examples/spectrogram_reconstruct.c` | Example: PNG -> Griffin-Lim -> WAV |
| `guides/spectrogram-reconstruct.md` | Guide with formulas + "Reading the formula in C" sections |

## Example App Flow

1. Read input WAV -> forward STFT -> write PNG spectrogram (round-trip demo)
2. Read PNG with `stbi_load()`
3. Invert Viridis colormap: nearest-neighbor in RGB space -> index `ci`
4. `dB = -80 + (ci / 255.0) * 80.0`, then `mag = pow(10, dB/20) * N`
5. Un-flip frequency axis (row 0 = highest freq in PNG)
6. `MD_griffin_lim(mag, ..., 50, signal_out, signal_len)`
7. `FIO_write_wav()` to output reconstructed audio

### Viridis Colormap Inversion
Brute-force nearest-neighbor in RGB space: for each pixel, compute `d = (r-lut_r)^2 + (g-lut_g)^2 + (b-lut_b)^2` against all 256 Viridis LUT entries. With 256 entries and typical image sizes, this is fast enough.

### Grayscale Detection
If `stbi_load` reports 1 channel, or R==G==B for all pixels, skip Viridis inversion and use gray value directly as colormap index.

## Key Tests

- **Round-trip**: signal -> STFT(mag+phase) -> iSTFT -> compare (tolerance 1e-6)
- **COLA property**: 50% overlap Hanning gives perfect reconstruction
- **Griffin-Lim convergence**: error decreases with more iterations
- **Frequency preservation**: reconstructed signal has energy at correct bins
- **Edge cases**: silence, single frame, single iteration

### Test Approach for Round-Trip
Since `MD_stft()` only returns magnitude, tests compute phase using public API:
```c
double window[N];
MD_Gen_Hann_Win(window, N);
for (unsigned f = 0; f < num_frames; f++) {
    double frame[N];
    for (unsigned n = 0; n < N; n++)
        frame[n] = signal[f * hop + n] * window[n];
    MD_magnitude_spectrum(frame, N, &mag[f * num_bins]);
    MD_phase_spectrum(frame, N, &phase[f * num_bins]);
}
```

## Verification

1. `make -C tests && ./tests/test_minidsp` — all tests pass
2. `make -C examples spectrogram_reconstruct` — builds cleanly
3. Run example on a generated spectrogram PNG, listen to output WAV
4. Compare forward STFT of reconstructed audio to original spectrogram visually

## Known Limitations

- **Griffin-Lim quality**: Produces audible metallic/phasey artifacts. Quality depends on STFT params matching, colormap inversion accuracy, and iteration count.
- **PNG quantization**: Viridis colormap reduces magnitude to 256 levels, introducing ~0.3 dB quantization noise.
- **Edge samples**: First and last partial-overlap regions have lower window normalization, producing slightly less accurate reconstruction at signal boundaries.
- **Memory**: Griffin-Lim needs ~`3 * num_frames * (N/2+1) * 8` bytes for working buffers. For 2s at 16kHz with N=1024, hop=16: ~24 MB.
