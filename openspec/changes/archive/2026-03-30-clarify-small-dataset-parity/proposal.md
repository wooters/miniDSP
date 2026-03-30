## Why

The VAD guide's comparison section mentions that miniDSP parameters were optimized on "train-clean-100" and the ViT was trained on "LibriSpeechConcat (small)", but it never explicitly states these are the same-sized dataset. This is an important detail: both systems used the small training split, making the comparison fair on data budget. Without this clarification, readers may assume the ViT had a data advantage.

## What Changes

- Update the "Systems under test" table in the comparison section to make it clear that both systems used the LibriVAD/LibriSpeechConcat **small** training data.
- Add a brief note near the table reinforcing that the data budgets are comparable, so the comparison isolates the modeling approach (hand-crafted features vs learned ViT) rather than data scale.

## Capabilities

### New Capabilities

_(none)_

### Modified Capabilities

- `vad-docs`: Updating the "Systems under test summary" requirement to specify that the guide must note both systems used the same-sized (small) training dataset.

## Impact

- **Docs**: Minor wording changes in `guides/vad.md` (comparison section only). No code or API changes.
