## Why

The mel_viz tool currently renders all frequency data using mel-scaled bands, but provides no way to see what a linear frequency scale looks like for comparison. Adding a mel/linear toggle would let users switch between the two in real time, making it immediately obvious how mel scaling better matches human perception — low frequencies get more visual resolution (where pitch discrimination is fine), while high frequencies are compressed (where perception is coarser). This turns the visualizer into a teaching tool, not just an aesthetic one.

## What Changes

- Add a "Scale" dropdown (or toggle) to the side panel with two options: **Mel** (default) and **Linear**
- In **Linear** mode, the frequency bands are uniformly spaced across 0–Nyquist rather than mel-spaced
- Both the **file mode** (C backend precomputation) and **mic mode** (JS live computation) support both scales
- The C backend generates both mel and linear energy arrays in `data.js` so switching is instant (no recomputation)
- The JS mic-mode provider builds both mel and linear filterbank weight matrices on initialization
- The renderer reads the current scale setting and uses the corresponding energy array each frame — no rendering changes needed beyond switching the data source

## Capabilities

### New Capabilities
- `linear-scale-mode`: Adds a linear (uniform Hz) frequency scale option alongside the existing mel scale, with a UI control to switch between them in real time

### Modified Capabilities

_(none — the mel pipeline is unchanged; linear mode is purely additive)_

## Impact

- **`tools/mel_viz/mel_viz.c`**: Compute and emit a second set of linear-spaced energy bands per frame
- **`tools/mel_viz/web/audio-provider.js`**: Add linear filterbank weights for mic mode; expose both mel and linear frame data
- **`tools/mel_viz/web/controls.js`**: Add Scale dropdown control
- **`tools/mel_viz/web/renderer.js`**: Read scale setting to select which energy array to render
- **`tools/mel_viz/web/index.html`**: Wire up the new control
- **miniDSP library (`src/`, `include/`)**: May need a new `MD_linear_energies()` function or the linear binning can be done directly in mel_viz.c using existing `MD_magnitude_spectrum()` / `MD_stft()` APIs
- **No breaking changes** — default behavior remains mel-scaled
