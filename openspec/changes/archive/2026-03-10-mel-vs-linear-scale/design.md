## Context

The mel_viz tool renders audio as radial concentric rings driven by mel-scaled frequency band energies. It operates in two modes: **file mode** (C backend precomputes mel energies into `data.js`, browser plays back) and **mic mode** (JS computes mel energies live via Web Audio API). The user wants a real-time toggle between mel-scaled and linear-scaled frequency bands so viewers can see why mel scaling better matches human perception.

Currently, both the C backend (`MD_mel_energies()`) and the JS mic provider (`buildMelWeights()`) produce mel-spaced band energies. The renderer (`renderer.js`) is scale-agnostic — it takes an array of band energies and maps them to rings. This means the renderer needs no changes beyond receiving the right data array.

## Goals / Non-Goals

**Goals:**
- Add a "Scale" dropdown to the side panel with "Mel" (default) and "Linear" options
- Switching scales is instant — no reloading, no recomputation delay
- Both file mode and mic mode support both scales
- The visual difference between mel and linear is immediately apparent (mel gives more resolution to low frequencies where human hearing is most discriminating)

**Non-Goals:**
- Adding a new `MD_linear_energies()` function to the miniDSP library — the linear binning is simple enough to do in mel_viz.c directly using `MD_magnitude_spectrum()`
- Changing the mel computation pipeline in any way
- Modifying the renderer's drawing logic
- Supporting arbitrary custom frequency scales (only mel and linear)

## Decisions

### 1. Linear energy computation strategy

**Decision**: Compute linear-spaced band energies by binning the one-sided PSD (from `MD_magnitude_spectrum()`) into uniformly-spaced frequency bands directly in `mel_viz.c`, rather than adding a library function.

**Rationale**: Linear binning is a straightforward summation of FFT bins into equal-width Hz bands. It doesn't warrant a library API — it's a ~15-line loop. The mel filterbank is a library function because the HTK mel mapping, triangular filter construction, and caching are non-trivial. Linear binning has no such complexity.

**Alternative considered**: Adding `MD_linear_energies()` to the library. Rejected because it would be a thin wrapper around magnitude spectrum with no reusable DSP value — it's just array slicing.

### 2. File mode: precompute both scale arrays in data.js

**Decision**: The C backend computes and emits both `frames` (mel energies) and `linearFrames` (linear energies) into `data.js`. The JS `FileProvider` exposes both and switches based on the current scale setting.

**Rationale**: Precomputing both avoids any runtime recomputation in the browser. The data.js file grows by ~2x in the frames section, but this is acceptable — for a typical 30s audio file at 30 fps with 24 bands, that's ~21K extra floats (~170 KB of text), well within browser memory limits.

**Alternative considered**: Computing linear energies in JS from the mel energies. Rejected — mel-to-linear is lossy (mel filters overlap and aren't invertible). The linear bands must be computed from the raw FFT in C.

### 3. Mic mode: dual filterbank weights

**Decision**: The JS `MicProvider` builds both mel and linear weight matrices at initialization. `getFrame()` checks the current scale setting and applies the corresponding weights.

**Rationale**: Building both weight matrices upfront is cheap (one-time cost). The per-frame cost is identical — same matrix-vector multiply, just different weights. The linear weights are trivial: each band spans `binCount / numBands` consecutive bins with weight 1.0 (rectangular, not triangular).

### 4. UI: dropdown control, not toggle button

**Decision**: Add a "Scale" dropdown to the controls panel (after "Rings") with options "Mel" and "Linear".

**Rationale**: A dropdown is consistent with the existing Palette and Rings controls. It's also extensible if other scales (e.g., Bark, ERB) are ever added, though that's a non-goal.

### 5. Settings integration

**Decision**: Add `scale: "mel"` to the `settings` object in `controls.js`. The `FileProvider` and `MicProvider` read `settings.scale` in their `getFrame()` methods to select which data to return.

**Rationale**: This follows the existing pattern — the renderer reads `settings` each frame, and the providers can do the same. The providers need access to `settings`, which means importing it from `controls.js`. The `MicProvider` already does its own computation per frame, so reading one extra field is negligible.

## Risks / Trade-offs

- **data.js size doubles** — The linear frames array adds ~170 KB for a 30s clip. This is fine for desktop browsers. If very long audio files become a concern, a future optimization could use binary format, but that's out of scope.
- **Bass envelope stays mel-derived** — The bass pulse is computed from the lowest mel bands regardless of scale setting. This is intentional: the bass pulse is a perceptual effect, and mel-derived bass tracking works well. Switching bass tracking to linear would make the pulse less responsive to perceived bass.
- **Linear scale will look "worse"** — By design. High-frequency bands dominate and low-frequency detail is lost, which is the whole pedagogical point.
