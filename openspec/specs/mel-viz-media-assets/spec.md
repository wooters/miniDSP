## Requirements

### Requirement: Compressed demo video committed to repo
The project SHALL include a compressed MP4 demo video of the mel_viz visualization at `docs/assets/mel-viz-demo.mp4`. The video SHALL be re-encoded from the source file to no more than 4 MB. The video SHALL include a compressed audio track (AAC) so the viewer can hear the music driving the visualization.

#### Scenario: Video file exists and is under size limit
- **WHEN** the file `docs/assets/mel-viz-demo.mp4` is checked
- **THEN** it SHALL exist and be no larger than 4 MB

#### Scenario: Video plays in browser
- **WHEN** the video is opened in Chrome, Firefox, or Safari
- **THEN** it SHALL play without errors (H.264 baseline profile, MP4 container)

### Requirement: Screenshot of mel_viz web UI committed to repo
The project SHALL include a PNG screenshot of the mel_viz visualization at `docs/assets/mel-viz-screenshot.png`. The screenshot SHALL show the radial visualization with visible rings and color.

#### Scenario: Screenshot file exists
- **WHEN** the file `docs/assets/mel-viz-screenshot.png` is checked
- **THEN** it SHALL exist and be a valid PNG image

### Requirement: Media assets registered in Doxyfile
Both `docs/assets/mel-viz-demo.mp4` and `docs/assets/mel-viz-screenshot.png` SHALL be listed in the Doxyfile `HTML_EXTRA_FILES` setting so Doxygen copies them to the output directory.

#### Scenario: Doxygen copies media to output
- **WHEN** `doxygen` is run
- **THEN** both `mel-viz-demo.mp4` and `mel-viz-screenshot.png` SHALL appear in the `docs/html/` output directory

### Requirement: Video embedded in guide page
The mel_viz guide page (`guides/mel-viz.md`) SHALL embed the demo video using `\htmlonly` with an HTML5 `<video>` tag. The video SHALL loop and show playback controls. It SHALL NOT autoplay (since it contains audio).

#### Scenario: Video visible in rendered Doxygen guide
- **WHEN** the mel_viz guide page is viewed in the Doxygen HTML output
- **THEN** the video SHALL be visible with playback controls, looping when played

### Requirement: Screenshot embedded in guide page
The mel_viz guide page SHALL embed the screenshot using Doxygen's `\image html` command.

#### Scenario: Screenshot visible in rendered Doxygen guide
- **WHEN** the mel_viz guide page is viewed in the Doxygen HTML output
- **THEN** the screenshot SHALL be displayed as a static image
