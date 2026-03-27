## Why

The `optimize/VAD/` directory contains a Python script (`optimize_vad.py`) that tunes VAD hyperparameters using Optuna. It depends on several third-party packages (`pyminidsp`, `optuna`, `numpy`, `scipy`, `soundfile`) but the miniDSP repo is a C library with no existing Python dependency management. There is currently no reproducible way to declare and install these dependencies.

## What Changes

- Add a `pyproject.toml` to `optimize/VAD/` declaring all Python dependencies with `requires-python >= 3.10`
- Generate and commit a `uv.lock` lockfile for reproducible installs
- Add `.venv/` and `__pycache__/` to the repo-root `.gitignore`
- Update `optimize/VAD/README.md` with `uv run` usage instructions

## Capabilities

### New Capabilities
- `vad-optimize-deps`: Self-contained uv project setup for VAD optimization script dependency management

### Modified Capabilities
<!-- No existing spec-level behavior changes -->

## Impact

- **New files**: `optimize/VAD/pyproject.toml`, `optimize/VAD/uv.lock`
- **Modified files**: `.gitignore`, `optimize/VAD/README.md`
- **Dependencies**: Introduces uv as a development tool prerequisite for running the optimization script
- **No impact** on the C library, its build system, or existing tests
