## Why

The VAD default parameters were hand-picked placeholders (equal 0.2 weights, threshold 0.5). A 300-trial Optuna optimization on LibriVAD (train-clean-100, all noise types, all SNRs) found parameters that raise F2 from **0.837 to 0.933** (+0.096). The optimized defaults should ship as the library's out-of-the-box configuration, and all documentation should explain how they were derived so users understand and trust the tuning.

## What Changes

- **Update `MD_vad_default_params()`** in `src/minidsp_vad.c` to use the optimized values from `optimize/VAD/best_params.json` (energy-dominant weighting, lower threshold, longer hangover, wider band).
- **Update the `vad-core` spec** to reflect the new default values as the normative reference.
- **Update Doxygen docs** in `include/minidsp.h` for `MD_vad_default_params()` to document the optimization provenance (dataset, metric, baseline vs. optimized scores).
- **Update the VAD tutorial guide** (`guides/vad.md`) with a new section explaining the optimization methodology, results, and per-condition breakdown.
- **Update the top-level README** (`README.md`) VAD optimization section to reference the results and link to the guide.
- **Update test expectations** in `tests/test_vad.c` if any tests assert specific default parameter values.

## Capabilities

### New Capabilities

_(none)_

### Modified Capabilities

- `vad-core`: Default parameter values change from hand-picked placeholders to Optuna-optimized values.
- `vad-docs`: Documentation updated to describe how defaults were derived (optimization methodology, dataset, metrics, results).

## Impact

- **`src/minidsp_vad.c`** — `MD_vad_default_params()` body changes.
- **`include/minidsp.h`** — Doxygen comment for `MD_vad_default_params()` updated.
- **`guides/vad.md`** — New "Default Parameter Optimization" section.
- **`README.md`** — VAD optimization section updated with results summary.
- **`tests/test_vad.c`** — Any assertions on default parameter values updated.
- **`openspec/specs/vad-core/spec.md`** — Default value requirements updated.
- **Behavioral change** — Users who rely on `MD_vad_default_params()` will get different (better) behavior after updating. This is intentional and not a breaking API change since the function signature and struct layout are unchanged.
