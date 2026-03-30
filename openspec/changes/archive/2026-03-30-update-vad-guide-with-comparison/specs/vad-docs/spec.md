## ADDED Requirements

### Requirement: Comparison section placement and structure
The guide `guides/vad.md` SHALL contain a `## Comparison with ViT-MFCC baseline` section placed after the "Default parameter optimization" section and before the "API summary" section. The section SHALL contain subsections for: systems under test, overall results, per-SNR breakdown, noise-type observations, dataset-specific optimization findings, interpretation guidance, and a pointer to the full comparison tooling.

#### Scenario: Section exists in correct position
- **WHEN** `guides/vad.md` is inspected
- **THEN** the "Comparison with ViT-MFCC baseline" heading SHALL appear after the "Re-tuning for your data" subsection and before the "API summary" heading

### Requirement: Systems under test summary
The comparison section SHALL include a summary table or description of the two systems under test: miniDSP VAD (hand-crafted features, C, 0 trainable parameters, 20ms frames) and ViT-MFCC small (Vision Transformer, ~3.5M parameters, 10ms frames). The summary SHALL note the evaluation dataset (LibriVAD test-clean, 702 files) and frame-rate handling (independent downsampling of labels to each system's native rate). The section SHALL include a link to the LibriVAD paper (https://arxiv.org/abs/2512.17281) so readers can find the full published results and methodology.

#### Scenario: Systems are clearly identified
- **WHEN** the comparison section is read
- **THEN** both systems SHALL be described with their type, parameter count, frame rate, and training/optimization data source

#### Scenario: Paper link is present
- **WHEN** the comparison section is read
- **THEN** it SHALL contain a hyperlink to the LibriVAD paper at https://arxiv.org/abs/2512.17281

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

### Requirement: Regenerate llms.txt
After `guides/vad.md` is updated, `llms.txt` and `llms-full.txt` SHALL be regenerated to include the new comparison content.

#### Scenario: llms.txt reflects comparison section
- **WHEN** `llms.txt` and `llms-full.txt` are inspected after regeneration
- **THEN** they SHALL include content from the new comparison section
