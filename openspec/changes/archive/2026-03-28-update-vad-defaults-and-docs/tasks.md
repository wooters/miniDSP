## 1. Update default parameters in code

- [x] 1.1 Update `MD_vad_default_params()` in `src/minidsp_vad.c` with optimized values from `best_params.json` (normalized weights, threshold, onset, hangover, adaptation rate, band frequencies)
- [x] 1.2 Update any test assertions in `tests/test_vad.c` that check specific default parameter values

## 2. Update vad-core spec

- [x] 2.1 Update `openspec/specs/vad-core/spec.md` — replace the default values in the "VAD params struct with sensible defaults" requirement with the optimized values table

## 3. Update Doxygen documentation

- [x] 3.1 Update `MD_vad_default_params()` doc-comment in `include/minidsp.h` — add `@note` block with optimization provenance (dataset, metric, trial count, baseline vs. optimized scores, pointer to guide)

## 4. Update tutorial guide

- [x] 4.1 Add "Default parameter optimization" section to `guides/vad.md` with: motivation, dataset description, methodology, baseline vs. optimized results table, key observations, per-condition summary, and re-tuning guidance

## 5. Update README

- [x] 5.1 Update the VAD optimization section in `README.md` with final results (F2=0.933, +0.096 improvement, dataset reference) and link to guide for full methodology

## 6. Verify

- [x] 6.1 Build the library (`make`) and run the test suite (`make test`) to confirm all tests pass with new defaults
- [x] 6.2 Regenerate `llms.txt` and `llms-full.txt` to include updated VAD content
