## Context

The VAD comparison script (`compare/VAD/compare_vad.py`) downloads a ViT-MFCC checkpoint from HuggingFace (`LibriVAD/LibriVAD` dataset repo) every time `download_vit_model()` is called. The `hf_hub_download` function does cache to the global `~/.cache/huggingface/hub/` directory, but this is invisible to the project and unreliable across environments. Users report repeated downloads on each comparison run.

## Goals / Non-Goals

**Goals:**
- Cache the ViT model checkpoint in a project-local directory so repeated runs skip the download.
- Keep model binaries out of git.

**Non-Goals:**
- Versioned model management or multi-model support.
- Overriding or disabling the global HuggingFace cache — the local cache simply adds a project-scoped layer.
- Offline-only mode or network error handling beyond what `hf_hub_download` already provides.

## Decisions

### 1. Cache location: `compare/VAD/.model_cache/`

Place the cache directory alongside the script that uses it, under a dot-prefixed name to signal it's generated.

**Alternatives considered:**
- **Top-level `.cache/` directory** — too generic; could conflict with other tools.
- **Rely on global HuggingFace cache** — current behavior; not portable across machines or CI.

### 2. Use `cache_dir` parameter of `hf_hub_download`

`hf_hub_download` already accepts a `cache_dir` kwarg that redirects its blob/symlink structure into a custom directory. This requires zero additional dependencies and preserves HuggingFace's built-in ETag-based freshness checking.

**Alternatives considered:**
- **Manual download + file existence check** — reinvents caching logic; loses ETag freshness.
- **`HUGGINGFACE_HUB_CACHE` env var** — global side-effect; would redirect all HF downloads, not just this model.

### 3. Single `.gitignore` entry at repo root

Add `compare/VAD/.model_cache/` to the repo-root `.gitignore` rather than adding a nested `.gitignore` inside `compare/VAD/`. Keeps ignore rules centralized.

## Risks / Trade-offs

- **Disk usage** — The cached model (~50 MB) persists locally. Acceptable for a dev/research project. Users can delete `.model_cache/` at any time to reclaim space.
- **Stale model** — If the upstream HuggingFace checkpoint is updated, `hf_hub_download` with `cache_dir` will still check ETags on each call (when online) and re-download if changed. No staleness risk in normal operation.
