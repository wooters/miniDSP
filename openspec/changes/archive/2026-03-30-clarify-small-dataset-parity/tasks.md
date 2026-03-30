## 1. Clarify small dataset parity in VAD guide

- [x] 1.1 Update the "Training data" row in the "Systems under test" table in `guides/vad.md` to label train-clean-100 as the small split (e.g., "train-clean-100 (small)")
- [x] 1.2 Add a sentence after the existing methodology prose (below the systems table) noting that both systems used the same-sized small training dataset, so the comparison isolates modeling approach rather than data scale

## 2. Regenerate llms.txt

- [x] 2.1 Regenerate `llms.txt` and `llms-full.txt` to include the updated wording
