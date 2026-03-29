## MODIFIED Requirements

### Requirement: Doxygen doc-comments for VAD API
All four public VAD functions (`MD_vad_default_params`, `MD_vad_init`, `MD_vad_calibrate`, `MD_vad_process_frame`) SHALL have full Doxygen doc-comments in `include/minidsp.h` with: LaTeX formula where applicable, `@param` for every parameter, `@return` description, `@code` usage example, and `@see` cross-references to `MD_energy`, `MD_zero_crossing_rate`, `MD_power_spectral_density`.

The `MD_vad_default_params` doc-comment SHALL additionally include a `@note` block stating: the defaults were optimized via 300-trial Optuna search on LibriVAD train-clean-100 (all noise types, all SNRs), optimizing F2 (beta=2); baseline F2=0.837 improved to F2=0.933 (P=0.782, R=0.981); and a pointer to `guides/vad.md` for full methodology.

#### Scenario: Doxygen generates complete API docs
- **WHEN** Doxygen is run on the project
- **THEN** all four VAD functions SHALL appear with complete documentation including formulas, parameter descriptions, return values, examples, and cross-references

#### Scenario: Default params doc includes optimization provenance
- **WHEN** the Doxygen output for `MD_vad_default_params` is viewed
- **THEN** it SHALL include a note describing the optimization dataset, metric, trial count, and baseline vs. optimized scores

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

## ADDED Requirements

### Requirement: README documents optimization results
The top-level `README.md` VAD optimization section SHALL include: the optimized F2 score (0.933), the improvement over baseline (+0.096), the dataset used (LibriVAD train-clean-100), and a link to the tutorial guide for full methodology.

#### Scenario: README VAD optimization section is current
- **WHEN** `README.md` is viewed
- **THEN** the VAD optimization section SHALL show the final optimization results (F2=0.933, baseline=0.837) and link to the guide
