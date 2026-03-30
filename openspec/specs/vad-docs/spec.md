## ADDED Requirements

### Requirement: Doxygen doc-comments for VAD API
All four public VAD functions (`MD_vad_default_params`, `MD_vad_init`, `MD_vad_calibrate`, `MD_vad_process_frame`) SHALL have full Doxygen doc-comments in `include/minidsp.h` with: LaTeX formula where applicable, `@param` for every parameter, `@return` description, `@code` usage example, and `@see` cross-references to `MD_energy`, `MD_zero_crossing_rate`, `MD_power_spectral_density`.

The `MD_vad_default_params` doc-comment SHALL additionally include a `@note` block stating: the defaults were optimized via 300-trial Optuna search on LibriVAD train-clean-100 (all noise types, all SNRs), optimizing F2 (beta=2); baseline F2=0.837 improved to F2=0.933 (P=0.782, R=0.981); and a pointer to `guides/vad.md` for full methodology.

#### Scenario: Doxygen generates complete API docs
- **WHEN** Doxygen is run on the project
- **THEN** all four VAD functions SHALL appear with complete documentation including formulas, parameter descriptions, return values, examples, and cross-references

#### Scenario: Default params doc includes optimization provenance
- **WHEN** the Doxygen output for `MD_vad_default_params` is viewed
- **THEN** it SHALL include a note describing the optimization dataset, metric, trial count, and baseline vs. optimized scores

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
The tutorial guide `guides/vad.md` SHALL contain: introduction, one `##` section per feature with LaTeX formula and "Reading the formula in C:" subsection, sections on adaptive normalization, weighted scoring, state machine (with transition table), interactive visualization iframe, and API summary with snippet embedding.

The guide SHALL additionally include a `## Default parameter optimization` section (placed after the state machine section and before the API summary) containing:
1. **Motivation** — why the defaults were optimized (placeholder values underperformed).
2. **Dataset** — LibriVAD train-clean-100, all noise types, all SNRs, 7560 files.
3. **Methodology** — Optuna TPE sampler, 300 trials, 10 parallel workers, F2 (beta=2) objective.
4. **Baseline vs. optimized** — table showing F2, precision, recall before and after.
5. **Key observations** — energy weight dominates, low threshold + long hangover favors recall, spectral entropy contributes little.
6. **Per-condition summary** — table or description of performance across noise types and SNR levels, noting worst-case (office noise at -5dB) and best-case.
7. **Guidance for re-tuning** — pointer to `optimize/VAD/` with instructions for running on custom data.

#### Scenario: Guide renders in Doxygen
- **WHEN** Doxygen is run on the project
- **THEN** the VAD guide page SHALL appear in the tutorials navigation with all formulas rendered, code snippets embedded, and visualization iframe functional

#### Scenario: Every formula has a C reading section
- **WHEN** the guide is inspected
- **THEN** every `\f[...\f]` display formula SHALL be followed by a "Reading the formula in C:" subsection with direct loop-level code

#### Scenario: Optimization methodology section is present
- **WHEN** the guide is inspected
- **THEN** it SHALL contain a "Default parameter optimization" section with dataset description, methodology, results table, key observations, per-condition summary, and re-tuning guidance

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

### Requirement: README documents optimization results
The top-level `README.md` VAD optimization section SHALL include: the optimized F2 score (0.933), the improvement over baseline (+0.096), the dataset used (LibriVAD train-clean-100), and a link to the tutorial guide for full methodology.

#### Scenario: README VAD optimization section is current
- **WHEN** `README.md` is viewed
- **THEN** the VAD optimization section SHALL show the final optimization results (F2=0.933, baseline=0.837) and link to the guide

### Requirement: Comparison section placement and structure
The guide `guides/vad.md` SHALL contain a `## Comparison with ViT-MFCC baseline` section placed after the "Default parameter optimization" section and before the "API summary" section. The section SHALL contain subsections for: systems under test, overall results, per-SNR breakdown, noise-type observations, dataset-specific optimization findings, interpretation guidance, and a pointer to the full comparison tooling.

#### Scenario: Section exists in correct position
- **WHEN** `guides/vad.md` is inspected
- **THEN** the "Comparison with ViT-MFCC baseline" heading SHALL appear after the "Re-tuning for your data" subsection and before the "API summary" heading

### Requirement: Systems under test summary
The comparison section SHALL include a summary table or description of the two systems under test: miniDSP VAD (hand-crafted features, C, 0 trainable parameters, 20ms frames) and ViT-MFCC small (Vision Transformer, ~3.5M parameters, 10ms frames). The summary SHALL note the evaluation dataset (LibriVAD test-clean, 702 files) and frame-rate handling (independent downsampling of labels to each system's native rate). The section SHALL include a link to the LibriVAD paper (https://arxiv.org/abs/2512.17281) so readers can find the full published results and methodology. The section SHALL explicitly note that both systems used the same-sized (small) training dataset, making the comparison fair on data budget and isolating the modeling approach as the key variable.

#### Scenario: Systems are clearly identified
- **WHEN** the comparison section is read
- **THEN** both systems SHALL be described with their type, parameter count, frame rate, and training/optimization data source

#### Scenario: Paper link is present
- **WHEN** the comparison section is read
- **THEN** it SHALL contain a hyperlink to the LibriVAD paper at https://arxiv.org/abs/2512.17281

#### Scenario: Small dataset parity is stated
- **WHEN** the comparison section is read
- **THEN** it SHALL explicitly state that both miniDSP (train-clean-100) and ViT-MFCC used the same-sized small training dataset

### Requirement: Overall results table
The comparison section SHALL include a results table showing F2, precision, recall, AUC (macro), AUC (pooled), and wall time for both miniDSP VAD and ViT-MFCC (small). The table SHALL use the default miniDSP parameters (train-clean-100 optimized).

#### Scenario: Overall metrics are present
- **WHEN** the overall results table is inspected
- **THEN** it SHALL show miniDSP F2=0.844, ViT F2=0.961, and the 22x speed difference

### Requirement: Per-SNR breakdown
The comparison section SHALL include a per-SNR breakdown table showing F2 and AUC for both systems across the six SNR levels (-5, 0, 5, 10, 15, 20 dB), with a note on how the F2 gap (~0.11-0.13) and AUC gap (~0.30-0.34) vary with SNR.

#### Scenario: SNR trends are shown
- **WHEN** the per-SNR table is inspected
- **THEN** it SHALL show both F2 and AUC values for each SNR level for both systems

### Requirement: Noise-type observations
The comparison section SHALL describe where miniDSP performs best (babble, SSN — distinctive spectral characteristics) and worst (transport, office — spectral overlap with speech) relative to the ViT baseline. It SHALL note the AUC vs F2 divergence: miniDSP's best F2 noise types (babble, SSN) have low AUC, indicating threshold-dependent performance.

#### Scenario: Best and worst noise types identified
- **WHEN** the noise-type discussion is read
- **THEN** it SHALL identify babble and SSN as miniDSP's strongest noise types and transport and office as the weakest, with explanation of why

### Requirement: Dataset-specific optimization findings
The comparison section SHALL summarize the dataset-specific optimization experiment: three parameter sets (train-clean-100, dev-clean, test-clean cheat) all achieve ~0.934-0.935 F2 on test-clean, confirming generalization. It SHALL state the F2 ceiling from threshold/weight tuning alone is ~0.935, and the remaining gap to ViT (0.935 → 0.961) requires richer features or a learned model.

#### Scenario: Generalization finding is stated
- **WHEN** the dataset-specific optimization subsection is read
- **THEN** it SHALL state that train-clean-100 parameters generalize to test-clean with negligible F2 loss and that the tuning ceiling is ~0.935

### Requirement: Interpretation guidance
The comparison section SHALL include guidance explaining: (1) AUC reveals score quality while F2 reflects a specific operating point, (2) miniDSP's low AUC (~0.65) vs ViT's high AUC (~0.97) means miniDSP's good F2 relies on Optuna threshold tuning rather than intrinsic score discriminability, and (3) practical implications — miniDSP is strong for embedded/real-time use with known conditions; ViT is superior when compute is available and conditions vary.

#### Scenario: AUC vs F2 interpretation is explained
- **WHEN** the interpretation guidance is read
- **THEN** it SHALL explain what the AUC gap means for score quality and threshold sensitivity

### Requirement: Reference to full comparison tooling
The comparison section SHALL include a pointer to `compare/VAD/` with a brief description of how to reproduce or extend the evaluation, and a reference to `compare/VAD/README_RESULTS.md` for complete per-condition breakdowns.

#### Scenario: Reproduction path is documented
- **WHEN** a user wants to reproduce the comparison
- **THEN** the section SHALL point them to `compare/VAD/` with enough context to get started

### Requirement: Regenerate llms.txt after comparison update
After `guides/vad.md` is updated, `llms.txt` and `llms-full.txt` SHALL be regenerated to include the new comparison content.

#### Scenario: llms.txt reflects comparison section
- **WHEN** `llms.txt` and `llms-full.txt` are inspected after regeneration
- **THEN** they SHALL include content from the new comparison section
