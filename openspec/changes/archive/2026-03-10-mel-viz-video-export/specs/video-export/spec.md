## ADDED Requirements

### Requirement: Export button visibility
The UI SHALL display an "Export" button in the controls area when the application is in file mode and the browser supports the WebCodecs API (`VideoEncoder` and `AudioEncoder`). The Export button SHALL NOT appear in mic mode or when the browser lacks WebCodecs support.

#### Scenario: File mode with supported browser (Chrome/Edge)
- **WHEN** the application loads in file mode in a browser that supports `VideoEncoder` and `AudioEncoder`
- **THEN** an "Export" button is visible in the controls area

#### Scenario: Mic mode
- **WHEN** the application is in mic mode
- **THEN** no Export button is displayed

#### Scenario: Unsupported browser
- **WHEN** the browser does not support `VideoEncoder` or `AudioEncoder`
- **THEN** no Export button is displayed

### Requirement: Browser compatibility note
The UI SHALL display a note near the Export button indicating that export requires Chrome or Edge.

#### Scenario: Export button visible
- **WHEN** the Export button is displayed
- **THEN** a text note reading "Export requires Chrome or Edge" (or similar) is visible near the button

### Requirement: Start export
When the user clicks the Export button, the system SHALL begin an offline export process that renders the visualization frame-by-frame from the precomputed mel data and encodes it with the decoded audio into an MP4 file. The export SHALL NOT require real-time audio playback. The Export button SHALL change to a "Cancel" button during export.

#### Scenario: User clicks Export
- **WHEN** the user clicks the Export button
- **THEN** the renderer's smoothing state is reset
- **AND** the system loops through all frames (0 to `numFrames - 1`), rendering each frame to the canvas and encoding it via `VideoEncoder`
- **AND** the audio file is decoded via `AudioContext.decodeAudioData()` and encoded via `AudioEncoder`
- **AND** video and audio are muxed into an MP4 container via mp4-muxer
- **AND** the Export button label changes to "Cancel"

### Requirement: Export completes with download
When the offline export finishes processing all frames, the system SHALL trigger a browser download of the resulting MP4 file.

#### Scenario: Export completes successfully
- **WHEN** all video frames and audio samples have been encoded and muxed
- **THEN** a download is triggered with filename `mel-viz-export.mp4`
- **AND** the button label reverts to "Export"

### Requirement: Cancel export
The user SHALL be able to cancel an in-progress export by clicking the Cancel button. Canceling SHALL discard the partial output and not trigger a download.

#### Scenario: User clicks Cancel during export
- **WHEN** the user clicks the Cancel button during an active export
- **THEN** the encoding and muxing are aborted
- **AND** no download is triggered
- **AND** the button label reverts to "Export"

### Requirement: Progress indication during export
The UI SHALL display a progress indicator during export showing how far along the export is (e.g., percentage or progress bar based on frames processed vs total frames).

#### Scenario: Export is active
- **WHEN** an export is in progress
- **THEN** a progress indicator (e.g., progress bar or percentage) is visible showing `framesProcessed / totalFrames`

### Requirement: Export uses MP4 format
The exported video SHALL be in MP4 container format with H.264 video and AAC audio codecs.

#### Scenario: Codec configuration
- **WHEN** the system initializes the `VideoEncoder`
- **THEN** it SHALL use the `avc1` (H.264) codec for video
- **AND** it SHALL use the `mp4a` (AAC) codec for audio via `AudioEncoder`
- **AND** video and audio SHALL be muxed into an MP4 container via mp4-muxer
