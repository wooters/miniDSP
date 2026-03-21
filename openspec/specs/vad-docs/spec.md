## ADDED Requirements

### Requirement: Doxygen doc-comments for VAD API
All four public VAD functions (`MD_vad_default_params`, `MD_vad_init`, `MD_vad_calibrate`, `MD_vad_process_frame`) SHALL have full Doxygen doc-comments in `include/minidsp.h` with: LaTeX formula where applicable, `@param` for every parameter, `@return` description, `@code` usage example, and `@see` cross-references to `MD_energy`, `MD_zero_crossing_rate`, `MD_power_spectral_density`.

#### Scenario: Doxygen generates complete API docs
- **WHEN** Doxygen is run on the project
- **THEN** all four VAD functions SHALL appear with complete documentation including formulas, parameter descriptions, return values, examples, and cross-references

### Requirement: VAD example program
An example program `examples/vad.c` SHALL demonstrate VAD usage: initialize with default params, optionally calibrate on silence frames, process frame-by-frame, output per-frame decision/score/features to CSV, and generate an interactive HTML visualization.

#### Scenario: Example builds and runs
- **WHEN** `make -C examples` is run
- **THEN** the `vad` example SHALL compile and link without errors

#### Scenario: Example produces CSV output
- **WHEN** the `vad` example is executed
- **THEN** it SHALL produce a CSV file with per-frame decision, score, and feature values

#### Scenario: Example produces HTML visualization
- **WHEN** the `vad` example is executed
- **THEN** it SHALL produce an interactive HTML file using Plotly.js

#### Scenario: Snippet markers present for guide embedding
- **WHEN** `examples/vad.c` is inspected
- **THEN** it SHALL contain `//! [vad-init]`, `//! [vad-calibrate]`, `//! [vad-process]`, and `//! [vad-custom-weights]` snippet marker pairs

### Requirement: Example Makefile integration
The `vad` example SHALL be added to `examples/Makefile` in both the `EXAMPLES` variable and the `plot:` target.

#### Scenario: Make plot runs vad
- **WHEN** `make -C examples plot` is run
- **THEN** the `vad` example SHALL be built and executed

### Requirement: Ignore file updates
`.gitignore` and `.dockerignore` SHALL each include entries for `examples/vad`, `examples/vad.csv`, and `examples/vad_plot.html`.

#### Scenario: Generated files are ignored
- **WHEN** the vad example is run from the examples directory
- **THEN** `git status` SHALL NOT show `vad`, `vad.csv`, or `vad_plot.html` as untracked files

### Requirement: Tutorial guide page
A tutorial guide `guides/vad.md` SHALL be created following existing guide conventions, containing: introduction, one `##` section per feature with LaTeX formula and "Reading the formula in C:" subsection, sections on adaptive normalization, weighted scoring, state machine (with transition table), interactive visualization iframe, and API summary with snippet embedding.

#### Scenario: Guide renders in Doxygen
- **WHEN** Doxygen is run on the project
- **THEN** the VAD guide page SHALL appear in the tutorials navigation with all formulas rendered, code snippets embedded, and visualization iframe functional

#### Scenario: Every formula has a C reading section
- **WHEN** the guide is inspected
- **THEN** every `\f[...\f]` display formula SHALL be followed by a "Reading the formula in C:" subsection with direct loop-level code

### Requirement: Interactive HTML visualization for docs
An interactive HTML visualization SHALL be created showing: waveform, per-frame feature traces, combined score with threshold line, and binary decision timeline. It SHALL use Plotly.js and target 380px iframe height.

#### Scenario: Visualization is iframe-embeddable
- **WHEN** the HTML file is embedded in an iframe in the Doxygen guide
- **THEN** it SHALL render correctly at 380px height without scrollbars

### Requirement: Guide navigation integration
The VAD guide SHALL be wired into documentation navigation: `guides/tutorials.md` SHALL include a `\subpage vad` entry, `include/minidsp.h` SHALL list VAD in the `@brief` feature list, and `README.md` SHALL mention VAD in the feature list.

#### Scenario: VAD appears in tutorials list
- **WHEN** the Doxygen tutorials page is viewed
- **THEN** a "Voice activity detection" subpage link SHALL appear

#### Scenario: README lists VAD
- **WHEN** `README.md` is viewed
- **THEN** VAD SHALL appear in the feature list

### Requirement: Doxyfile updates for VAD assets
The `Doxyfile` SHALL be updated to include any generated HTML visualization files in `HTML_EXTRA_FILES`.

#### Scenario: Assets are copied to docs output
- **WHEN** Doxygen is run
- **THEN** all VAD-related HTML assets SHALL be available in the docs output directory

### Requirement: Regenerate llms.txt
After all code and documentation changes are complete, `llms.txt` and `llms-full.txt` SHALL be regenerated to include the new VAD content.

#### Scenario: llms.txt includes VAD
- **WHEN** `llms.txt` and `llms-full.txt` are inspected after regeneration
- **THEN** they SHALL reference the new VAD API and guide content
