## Context

The comparison section's "Systems under test" table lists training data as "Optuna optimization on train-clean-100" vs "Supervised training on LibriSpeechConcat (small)". These are the same data split (small), but that equivalence isn't stated. The table's "Training data" row and the surrounding prose both need a small tweak.

## Goals / Non-Goals

**Goals:**
- Make it explicit that both systems used the same-sized (small) training data, so the comparison is fair on data budget.

**Non-Goals:**
- No changes to the optimization section, parameter generalization section, or any other part of the guide.
- No changes to the comparison results or methodology.

## Decisions

### 1. Annotate "train-clean-100" as the small split in the table

Update the Training data row to read "Optuna optimization on train-clean-100 (small)" to mirror the ViT's "LibriSpeechConcat (small)" label. This makes the parity visible at a glance.

### 2. Add a one-sentence note below the table

After the existing prose about frame-rate handling and thresholds, add a sentence stating that both systems used the same-sized training data, isolating the modeling approach as the key variable.

## Risks / Trade-offs

None — this is a minor wording clarification with no risk.
