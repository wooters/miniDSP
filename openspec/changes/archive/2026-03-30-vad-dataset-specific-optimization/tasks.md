## 1. Run Optimizations

- [x] 1.1 Run Optuna optimization on `dev-clean` split (300 trials), save results to `optimize/VAD/best_params_dev_clean.json`
- [x] 1.2 Run Optuna optimization on `test-clean` split (300 trials), save results to `optimize/VAD/best_params_test_clean.json`

## 2. Extend Comparison Script

- [x] 2.1 Add `load_vad_params(json_path)` helper function that reads a `best_params.json` file and returns a kwargs dict for `VAD()`
- [x] 2.2 Add `--params` repeatable CLI argument in `label=path` format, with fallback to library defaults when omitted
- [x] 2.3 Modify `eval_minidsp()` to accept optional `vad_kwargs` dict and pass to `VAD()` constructor
- [x] 2.4 Update `main()` to loop over all parameter sets, running `eval_minidsp()` for each
- [x] 2.5 Update `print_overall()` to display one row per miniDSP configuration
- [x] 2.6 Update `print_breakdown()` and `_print_bucket_table()` to display one column per miniDSP configuration

## 3. Verify

- [x] 3.1 Run comparison with all three parameter sets and verify output tables display correctly
- [x] 3.2 Verify backward compatibility: running without `--params` produces the same output as before
