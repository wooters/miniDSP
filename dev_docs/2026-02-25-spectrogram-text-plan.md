# Spectrogram Text Art — Implementation Plan

## Task 1: Bitmap font data

Create `src/minidsp_spectext.c` with a static 5x7 bitmap font covering printable ASCII (32–126). Each character is stored as 5 bytes (one per column), each byte encoding 7 row bits. Include a helper function that, given a string, returns the pixel state at any (column, row) coordinate.

**Files:** `src/minidsp_spectext.c`

---

## Task 2: Core synthesis function

Implement `MD_spectrogram_text()` in `src/minidsp_spectext.c`:

1. Compute bitmap grid dimensions: `W = len * 6 - 1` columns, 7 rows.
2. Compute `col_samples = (duration_sec / W) * sample_rate`.
3. For each column, sum sine waves for all "on" pixels. Carry phase across columns for continuous tones.
4. Apply 3 ms raised-cosine fade on tone onsets/offsets between adjacent columns.
5. Normalize output to 0.9 peak amplitude.
6. Return number of samples written.

Assert preconditions: non-null pointers, `freq_lo < freq_hi`, `freq_hi <= sample_rate / 2`, `duration_sec > 0`, `sample_rate > 0`, non-empty text, `max_len` sufficient.

**Files:** `src/minidsp_spectext.c`

---

## Task 3: Public header and build integration

- Add `MD_spectrogram_text()` declaration and Doxygen doc-comment to `include/minidsp.h`.
- Add "spectrogram text art" to the `@brief` feature list in `minidsp.h`.
- Add `src/minidsp_spectext.o` to the root `Makefile` object list.

**Files:** `include/minidsp.h`, `Makefile`

---

## Task 4: Tests

Add test cases to `tests/test_minidsp.c`:

- Output length matches expected `col_samples * W`.
- Single character "A": verify non-zero energy in the output.
- Frequency check: for a single column with one "on" pixel, run `MD_magnitude_spectrum` and verify the peak is at the expected frequency bin.
- Normalization: verify peak absolute value is approximately 0.9.
- Space character produces silence (all zeros).

**Files:** `tests/test_minidsp.c`

---

## Task 5: Fetch stb_image_write.h

Download `stb_image_write.h` (public domain, single-header PNG writer) into a `third_party/` directory. This header is used only by the example program, not the library.

**Files:** `third_party/stb_image_write.h`

---

## Task 6: Example program

Create `examples/spectrogram_text.c`:

1. Parse command-line arguments: text string, optional `--colormap viridis|grayscale` (default: viridis).
2. Call `MD_spectrogram_text()` with defaults: `freq_lo=200`, `freq_hi=8000`, `duration=2.0`, `sample_rate=16000`.
3. Write WAV file via libsndfile.
4. Compute STFT with large FFT (e.g., 1024) and small hop (e.g., 16) for a sharp spectrogram.
5. Write HTML file with embedded Plotly spectrogram.
6. Write PNG file using `stb_image_write.h` with the selected colormap (Viridis LUT or grayscale).

Viridis colormap: embed a 256-entry RGB lookup table. Grayscale: linear map from dB magnitude to 0–255.

**Files:** `examples/spectrogram_text.c`, `examples/Makefile`

---

## Task 7: Build and ignore file updates

- Add example to `EXAMPLES` and `plot` target in `examples/Makefile`.
- Add three lines to `.gitignore` and `.dockerignore`: `examples/spectrogram_text`, `examples/spectrogram_text.csv`, `examples/spectrogram_text.html`.
- Also add `examples/spectrogram_text.png` and `examples/spectrogram_text.wav` to both ignore files.

**Files:** `examples/Makefile`, `.gitignore`, `.dockerignore`

---

## Task 8: Guide page

Create `guides/spectrogram-text.md`:

- Explain the concept: bitmap font → frequency mapping → sine synthesis.
- Show the formula mapping with a "Reading the formula in C" section.
- Embed the HTML spectrogram via iframe.
- Show the generated PNG image.
- Add `\subpage spectrogram-text` to `guides/tutorials.md`.

**Files:** `guides/spectrogram-text.md`, `guides/tutorials.md`

---

## Task 9: Doxyfile and README

- Add generated HTML assets to `HTML_EXTRA_FILES` in `Doxyfile`.
- Add spectrogram text to the feature list in `README.md`.

**Files:** `Doxyfile`, `README.md`

---

## Dependency Order

Tasks 1–2 are sequential (font before synthesis). Task 3 depends on Task 2 (need the function to declare). Task 4 depends on Task 3 (tests need the header). Task 5 is independent. Task 6 depends on Tasks 3 and 5. Task 7 depends on Task 6. Task 8 depends on Task 6. Task 9 depends on Task 8.

```
1 → 2 → 3 → 4
              ↘
    5 --------→ 6 → 7
                 ↘
                  8 → 9
```
