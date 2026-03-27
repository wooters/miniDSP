## 1. Project Configuration

- [x] 1.1 Create `optimize/VAD/pyproject.toml` with project name `minidsp-optimize-vad`, version `0.1.0`, `requires-python = ">=3.10"`, and dependencies: `pyminidsp`, `optuna`, `soundfile`, `numpy`, `scipy`
- [x] 1.2 Run `uv lock` in `optimize/VAD/` to generate `uv.lock`

## 2. Git Hygiene

- [x] 2.1 Append `.venv/` and `__pycache__/` entries to the repo-root `.gitignore`

## 3. Documentation

- [x] 3.1 Update `optimize/VAD/README.md` to replace bare `python` invocations with `uv run python` equivalents and add a prerequisites note about installing uv

## 4. Verification

- [x] 4.1 Run `uv run python optimize_vad.py --help` from `optimize/VAD/` and confirm it resolves dependencies and prints help without errors
