## Why

The current VAD comparison uses miniDSP parameters optimized on `train-clean-100`, which is the standard practice (optimize on training data, evaluate on test data). However, we want to understand how much performance varies by optimization dataset. Adding `dev-clean`-optimized and `test-clean`-optimized ("cheating") parameter sets to the comparison will reveal: (1) whether the train-clean-100 parameters generalize well or overfit to that split, and (2) the upper bound of what threshold/weight tuning alone can achieve when perfectly matched to the evaluation data.

## What Changes

- Run Optuna optimization on `dev-clean` split, producing a second set of optimized VAD parameters
- Run Optuna optimization on `test-clean` split (the same data used for evaluation), producing a "cheating" parameter set that represents the best possible tuning
- Extend `compare_vad.py` to evaluate multiple miniDSP VAD configurations side by side (train-clean-100, dev-clean, test-clean-optimized) against the ViT-MFCC baseline
- Store each parameter set as a separate JSON file with clear naming

## Capabilities

### New Capabilities

- `multi-param-comparison`: Extend the VAD comparison framework to evaluate multiple miniDSP parameter sets in a single run, with labeled columns in the output tables

### Modified Capabilities

## Impact

- `compare/VAD/compare_vad.py` — Must support multiple miniDSP parameter configurations instead of just the library defaults
- `optimize/VAD/` — New optimization runs producing `best_params_dev_clean.json` and `best_params_test_clean.json`
- No changes to the C library itself; all parameter sets are applied at the Python level via `pyminidsp.VAD()` constructor
- No changes to the ViT-MFCC baseline evaluation
