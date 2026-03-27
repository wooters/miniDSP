# VAD Optimization uv Setup -- Implementation Plan

## Tasks

### 1. Create `optimize/VAD/pyproject.toml`

```toml
[project]
name = "minidsp-optimize-vad"
version = "0.1.0"
description = "Hyperparameter optimization for miniDSP VAD"
requires-python = ">=3.10"
dependencies = [
    "pyminidsp",
    "optuna",
    "soundfile",
    "numpy",
    "scipy",
]
```

No `[build-system]` section -- this is an application, not an installable
package.

### 2. Generate lockfile

Run `uv lock` inside `optimize/VAD/` to produce `uv.lock`. Commit the
lockfile so builds are reproducible.

### 3. Add `.venv/` to repo-root `.gitignore`

Append a Python section with `.venv/` so uv's auto-created virtual
environments are never committed. Also add `__pycache__/`.

### 4. Update `optimize/VAD/README.md`

Replace bare `python` invocations with `uv run python` equivalents. Add a
brief prerequisites note mentioning uv installation.

### 5. Verify

Run `uv run python optimize_vad.py --help` from `optimize/VAD/` to confirm
dependency resolution and script execution work.
