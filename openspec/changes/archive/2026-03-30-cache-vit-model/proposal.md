## Why

Every run of the VAD comparison script (`compare/VAD/compare_vad.py`) calls `hf_hub_download()` without a local `cache_dir`, relying on the user's global HuggingFace cache. This means the model is re-downloaded on fresh environments, CI, or whenever the global cache is cleared. A project-local cache makes runs reproducible and avoids repeated network hits.

## What Changes

- Add a local model cache directory under `compare/VAD/` (e.g., `compare/VAD/.model_cache/`) where downloaded HuggingFace checkpoints are stored.
- Update `download_vit_model()` to pass `cache_dir` pointing to the local cache directory, so `hf_hub_download` reuses previously downloaded files.
- Add the cache directory to `.gitignore` so model binaries are never committed.

## Capabilities

### New Capabilities
- `model-cache`: Local caching of the ViT-MFCC HuggingFace checkpoint so it persists across runs without re-downloading.

### Modified Capabilities
_(none — no existing spec-level requirements change)_

## Impact

- **Code**: `compare/VAD/compare_vad.py` — `download_vit_model()` function gains a `cache_dir` parameter.
- **Files**: New `.model_cache/` directory created at runtime under `compare/VAD/`.
- **Git**: `.gitignore` updated with `compare/VAD/.model_cache/` entry.
- **Dependencies**: No new dependencies — `hf_hub_download` already supports `cache_dir`.
