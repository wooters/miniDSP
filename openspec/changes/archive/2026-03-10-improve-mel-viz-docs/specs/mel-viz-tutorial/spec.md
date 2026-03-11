## ADDED Requirements

### Requirement: Step-by-step tutorial section in guide
The mel_viz guide page SHALL include a "Visualize your own audio" section with numbered steps that walk the reader from a WAV file to a running browser visualization.

#### Scenario: Tutorial covers build step
- **WHEN** the reader follows the tutorial
- **THEN** the first step SHALL instruct them to build mel_viz from the repository root using `make tools`

#### Scenario: Tutorial covers running mel_viz on a WAV file
- **WHEN** the reader has built mel_viz
- **THEN** a step SHALL show the command to run mel_viz on their own WAV file with an output directory flag

#### Scenario: Tutorial covers serving the output
- **WHEN** the reader has generated the output folder
- **THEN** a step SHALL show how to start a local HTTP server and open the visualization in a browser

#### Scenario: Tutorial covers customizing visual controls
- **WHEN** the reader has the visualization running in a browser
- **THEN** a step SHALL describe the available visual controls (palette, smoothing, bass, wobble, glow, rings) and how to adjust them in the side panel

#### Scenario: Tutorial covers video export
- **WHEN** the reader wants to share their visualization
- **THEN** a step SHALL describe how to export an MP4 video using the built-in export feature (Chrome/Edge only)

### Requirement: Guide page leads with visuals
The guide page SHALL be restructured so that the video demo and screenshot appear before the tutorial and reference sections. The reader SHALL see what mel_viz produces before being asked to build or run anything.

#### Scenario: Visual content precedes tutorial
- **WHEN** the guide page is read top to bottom
- **THEN** the video and screenshot SHALL appear before the "Visualize your own audio" tutorial section

### Requirement: Guide page ends with reference material
Architecture diagrams, CLI flags, and technical details SHALL appear after the tutorial section, not before it.

#### Scenario: Architecture section is at the end
- **WHEN** the guide page is read top to bottom
- **THEN** the architecture diagram and control reference SHALL appear after the tutorial
