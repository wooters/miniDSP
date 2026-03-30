## Context

The VAD comparison framework (`compare/VAD/compare_vad.py`) currently evaluates a single miniDSP VAD configuration (library defaults, which are the train-clean-100-optimized parameters) against the ViT-MFCC baseline. The optimization script (`optimize/VAD/optimize_vad.py`) already supports arbitrary dataset/split selection via `--split`. The `pyminidsp.VAD()` constructor accepts all tuning parameters as keyword arguments, and `best_params.json` stores the full parameter set including `weights_normalized`.

## Goals / Non-Goals

**Goals:**
- Run optimization on `dev-clean` and `test-clean` splits, producing two new parameter JSON files
- Extend `compare_vad.py` to load and evaluate multiple miniDSP parameter sets in a single run
- Display all miniDSP configurations alongside ViT-MFCC in the output tables
- Keep the same evaluation methodology (same test-clean data, same metrics, same beta)

**Non-Goals:**
- Changing the C library defaults (they remain train-clean-100-optimized)
- Modifying the optimization script itself (it already supports `--split`)
- Adding new metrics or changing the evaluation methodology
- Optimizing the ViT-MFCC baseline

## Decisions

### 1. Parameter loading: JSON files with a naming convention

Each optimization run produces a `best_params_<split>.json` file in `optimize/VAD/`. The comparison script loads them by path.

**Rationale:** The optimizer already outputs JSON in a consistent format. Using separate files with a clear naming convention (e.g., `best_params_dev_clean.json`, `best_params_test_clean.json`) keeps things simple and avoids modifying the optimizer. The existing `best_params.json` is the train-clean-100 result.

**Alternative considered:** A single JSON file with all parameter sets keyed by split name. Rejected because each optimization run is independent and the optimizer writes one file at a time.

### 2. Comparison script: `--params` flag to specify named parameter sets

Add a repeatable `--params` CLI argument to `compare_vad.py`:
```
--params "train-clean-100=../optimize/VAD/best_params.json" \
--params "dev-clean=../optimize/VAD/best_params_dev_clean.json" \
--params "test-clean (cheat)=../optimize/VAD/best_params_test_clean.json"
```

Each entry is `label=path`. When no `--params` is provided, fall back to the current behavior (library defaults, labeled "miniDSP VAD").

**Rationale:** This is flexible — any number of parameter sets, any labels, any JSON files. It doesn't hardcode split names or file paths. The label appears in output tables.

**Alternative considered:** Auto-discovering JSON files by glob pattern. Rejected because explicit is better — the user controls which configurations appear and what they're called.

### 3. VAD instantiation: load JSON and pass to `VAD()` constructor

Create a helper function `load_vad_params(json_path) -> dict` that reads a `best_params.json` file and returns the kwargs dict for `VAD()`:
```python
{
    "weights": normalized_weights,
    "threshold": ...,
    "onset_frames": ...,
    "hangover_frames": ...,
    "adaptation_rate": ...,
    "band_low_hz": ...,
    "band_high_hz": ...,
}
```

The existing `eval_minidsp()` function gets a new optional `vad_kwargs` parameter. When provided, it passes them to `VAD(**vad_kwargs)` instead of using bare `VAD()`.

**Rationale:** Minimal change to existing code. The JSON format is already established by the optimizer.

### 4. Output format: extend existing tables with additional columns

The overall table and breakdown tables add one column per miniDSP parameter set. The ViT-MFCC column remains as-is. Each miniDSP configuration gets its own labeled column.

**Rationale:** Side-by-side comparison in a single table is the most readable format for this use case.

## Risks / Trade-offs

- **Runtime scales linearly** with number of parameter sets — each miniDSP configuration does a full pass over all test files. Mitigation: miniDSP evaluation is fast (~3.6s per run), so even 3 configurations is under 15s total.
- **JSON format coupling** — the comparison script depends on the optimizer's output format. Mitigation: the format is simple and stable; the helper function isolates the parsing.
- **test-clean optimization is "cheating"** — results will be optimistic by design. Mitigation: label it clearly in output (e.g., "test-clean (cheat)") so it's never mistaken for a legitimate result.
