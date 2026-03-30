## MODIFIED Requirements

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
