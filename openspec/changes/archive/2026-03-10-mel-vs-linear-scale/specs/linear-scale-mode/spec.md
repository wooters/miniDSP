## ADDED Requirements

### Requirement: Scale dropdown control
The mel_viz side panel SHALL include a "Scale" dropdown with options "Mel" (default) and "Linear", placed after the "Rings" dropdown.

#### Scenario: Default scale is mel
- **WHEN** the visualizer loads
- **THEN** the Scale dropdown SHALL display "Mel" and the visualization SHALL use mel-scaled frequency bands

#### Scenario: User switches to linear
- **WHEN** the user selects "Linear" from the Scale dropdown
- **THEN** the visualization SHALL immediately switch to linear-scaled frequency bands without interrupting audio playback

#### Scenario: User switches back to mel
- **WHEN** the user selects "Mel" from the Scale dropdown after viewing linear
- **THEN** the visualization SHALL immediately switch back to mel-scaled frequency bands

### Requirement: Linear frequency band computation in file mode
The C backend SHALL compute linear-spaced frequency band energies alongside mel energies and include both in the output `data.js` file.

#### Scenario: data.js contains both frame arrays
- **WHEN** the C backend processes a WAV file
- **THEN** the output `data.js` SHALL contain both a `frames` array (mel energies) and a `linearFrames` array (linear energies), each with `numFrames * numBands` values

#### Scenario: Linear bands are uniformly spaced
- **WHEN** computing linear band energies for a frame
- **THEN** the frequency range `[min_freq, max_freq]` SHALL be divided into `numBands` equal-width Hz bands, and each band's energy SHALL be the sum of PSD values for FFT bins falling within that band's frequency range

#### Scenario: Linear and mel arrays have identical dimensions
- **WHEN** the backend generates `data.js`
- **THEN** `frames` and `linearFrames` SHALL have the same length (`numFrames * numBands`) so the frontend can index them identically

### Requirement: Linear frequency band computation in mic mode
The JS `MicProvider` SHALL support both mel and linear filterbank weights for live microphone input.

#### Scenario: Mic mode builds both weight matrices
- **WHEN** `MicProvider.start()` initializes
- **THEN** it SHALL build both mel filterbank weights and linear filterbank weights for the AnalyserNode's bin count

#### Scenario: Mic mode switches scales in real time
- **WHEN** the Scale setting changes during live mic input
- **THEN** `getFrame()` SHALL return energies computed with the weight matrix corresponding to the current scale setting

### Requirement: File mode provider supports both scales
The JS `FileProvider` SHALL return either mel or linear frame data based on the current scale setting.

#### Scenario: FileProvider returns mel data by default
- **WHEN** `getFrame()` is called with Scale set to "Mel"
- **THEN** it SHALL return the mel energy slice from `data.frames`

#### Scenario: FileProvider returns linear data when selected
- **WHEN** `getFrame()` is called with Scale set to "Linear"
- **THEN** it SHALL return the linear energy slice from `data.linearFrames`

### Requirement: Bass envelope is scale-independent
The bass envelope SHALL always be derived from mel-scaled bands, regardless of the current scale setting.

#### Scenario: Bass pulse unchanged in linear mode
- **WHEN** the Scale dropdown is set to "Linear"
- **THEN** the bass envelope (`getBass()`) SHALL still be computed from the lowest mel bands, not from linear bands
