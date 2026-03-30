## 1. Git Configuration

- [x] 1.1 Add `compare/VAD/.model_cache/` entry to the repo-root `.gitignore`

## 2. Update Download Function

- [x] 2.1 Define `CACHE_DIR` constant (relative to the script's directory) pointing to `.model_cache/`
- [x] 2.2 Update `download_vit_model()` to pass `cache_dir=CACHE_DIR` to `hf_hub_download`

## 3. Verification

- [x] 3.1 Run `compare_vad.py` once and confirm the model is saved under `compare/VAD/.model_cache/`
- [x] 3.2 Run again and confirm no re-download occurs (check output or network activity)
- [x] 3.3 Confirm `git status` does not show `.model_cache/` files
