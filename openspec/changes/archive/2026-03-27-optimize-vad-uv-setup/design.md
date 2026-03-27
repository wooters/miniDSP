## Context

The miniDSP repo is a C library. The `optimize/VAD/` directory contains `optimize_vad.py`, a Python script that tunes VAD hyperparameters via Optuna. It depends on `pyminidsp`, `optuna`, `numpy`, `scipy`, and `soundfile`. There is no Python dependency management in the repo today -- users must manually install packages into whatever Python environment they happen to have active.

## Goals / Non-Goals

**Goals:**
- Declare Python dependencies in a standard `pyproject.toml` so they are explicit and versioned
- Pin exact dependency versions via `uv.lock` for reproducible environments
- Enable single-command setup and execution via `uv run`
- Keep the Python tooling fully isolated within `optimize/VAD/`

**Non-Goals:**
- Adding Python dependency management at the repo root
- Packaging the optimization script as a distributable Python package
- Pinning or vendoring uv itself
- Migrating other scripts or tools to uv

## Decisions

### 1. Self-contained uv project in `optimize/VAD/`

**Choice**: Place `pyproject.toml` and `uv.lock` directly in `optimize/VAD/` rather than at the repo root or in `optimize/`.

**Rationale**: Each optimization tool may have different dependencies. Co-locating the project file with the script avoids cross-tool dependency conflicts and keeps the C project root clean.

**Alternatives considered**:
- Repo-root `pyproject.toml`: Would mix Python metadata into a C project and create a single dependency set for all future Python tools.
- `requirements.txt`: No standard lockfile, no automatic venv management, less reproducible.
- Docker/devcontainer: Higher overhead for a single script; uv is lighter-weight.

### 2. Application-only `pyproject.toml` (no `[build-system]`)

**Choice**: Omit `[build-system]` since this is a runnable script, not an installable package.

**Rationale**: Adding build metadata would imply the script is pip-installable, which it is not. uv handles application-style projects without a build backend.

### 3. Commit `uv.lock`

**Choice**: Commit the lockfile to version control.

**Rationale**: Ensures all developers and CI get identical dependency versions. This is standard practice for application-style projects.

### 4. `.venv/` in repo-root `.gitignore`

**Choice**: Add `.venv/` to the repo-root `.gitignore` rather than a local `optimize/VAD/.gitignore`.

**Rationale**: Any future Python tool under `optimize/` will also produce a `.venv/`. A single root-level ignore covers all of them.

## Risks / Trade-offs

- **[uv not installed]** -> Users must install uv separately. Mitigated by documenting installation in the README.
- **[Lock staleness]** -> `uv.lock` can drift from `pyproject.toml` if someone adds a dependency but forgets to re-lock. Mitigated by `uv run` auto-checking consistency.
- **[pyminidsp compatibility]** -> `pyminidsp` must be installable from PyPI (or a configured index). If it requires a local build of the C library, that dependency chain needs documentation. Mitigated by verifying with `uv run python optimize_vad.py --help`.
