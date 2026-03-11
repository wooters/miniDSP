## 1. C Backend — Linear Energy Computation

- [x] 1.1 Add `compute_linear_frames()` function to `mel_viz.c` that bins `MD_magnitude_spectrum()` output into `numBands` uniformly-spaced Hz bands over `[min_freq, max_freq]`
- [x] 1.2 Call `compute_linear_frames()` in `main()` alongside `compute_mel_frames()` and pass both arrays to `assemble_output()`
- [x] 1.3 Update `write_data_js()` to emit the `linearFrames` array in `data.js` alongside the existing `frames` array

## 2. JS Frontend — FileProvider

- [x] 2.1 Update `FileProvider` constructor to store `data.linearFrames`
- [x] 2.2 Update `FileProvider.getFrame()` to read `settings.scale` and return either mel or linear slice

## 3. JS Frontend — MicProvider

- [x] 3.1 Add `buildLinearWeights()` function to `audio-provider.js` that creates a rectangular filterbank with uniformly-spaced Hz bands
- [x] 3.2 Update `MicProvider.start()` to build both mel and linear weight matrices
- [x] 3.3 Update `MicProvider.getFrame()` to select weight matrix based on `settings.scale`

## 4. UI Controls

- [x] 4.1 Add `scale: "mel"` to the `settings` object in `controls.js`
- [x] 4.2 Add a "Scale" dropdown (options: "Mel", "Linear") to `initControls()`, placed after the "Rings" dropdown

## 5. Integration and Testing

- [x] 5.1 Rebuild mel_viz, regenerate output from a test WAV, and verify `data.js` contains both `frames` and `linearFrames` with matching dimensions
- [x] 5.2 Open in browser and verify the Scale dropdown toggles between mel and linear visualization in file mode
- [x] 5.3 Test mic mode: verify both scales work with live microphone input
- [x] 5.4 Verify bass envelope remains mel-derived regardless of scale setting
