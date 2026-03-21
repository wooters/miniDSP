## ADDED Requirements

### Requirement: VAD visualization filename must not collide with Doxygen page output
The VAD example program SHALL write its HTML visualization to `vad_plot.html` (not `vad.html`) to avoid overwriting the Doxygen-generated guide page in the docs output directory.

#### Scenario: No filename collision in docs build
- **WHEN** Doxygen builds the documentation and copies `HTML_EXTRA_FILES` to `docs/html/`
- **THEN** the guide page (`vad.html`, generated from `guides/vad.md` with anchor `{#vad}`) and the visualization (`vad_plot.html`, from `HTML_EXTRA_FILES`) SHALL both exist as separate files

#### Scenario: Iframe displays the visualization
- **WHEN** a user views the VAD guide page in the rendered documentation
- **THEN** the iframe in the "Visualization" section SHALL load `vad_plot.html` and display the interactive Plotly chart (not a recursive embed of the guide page)
