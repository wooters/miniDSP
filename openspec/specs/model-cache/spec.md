### Requirement: Local model cache directory
The system SHALL store downloaded HuggingFace model checkpoints in a project-local cache directory at `compare/VAD/.model_cache/`.

#### Scenario: Cache directory created on first download
- **WHEN** `download_vit_model()` is called and `.model_cache/` does not exist
- **THEN** the directory SHALL be created automatically before the download proceeds

#### Scenario: Cache directory already exists
- **WHEN** `download_vit_model()` is called and `.model_cache/` already exists
- **THEN** the function SHALL reuse the existing directory without error

### Requirement: Download uses local cache
The `download_vit_model()` function SHALL pass `cache_dir` pointing to `compare/VAD/.model_cache/` when calling `hf_hub_download`, so that the downloaded checkpoint is stored locally and reused on subsequent runs.

#### Scenario: First download populates local cache
- **WHEN** `download_vit_model()` is called for the first time (empty cache)
- **THEN** the model file SHALL be downloaded from HuggingFace and stored under `.model_cache/`

#### Scenario: Subsequent run uses cached model
- **WHEN** `download_vit_model()` is called and the model already exists in `.model_cache/`
- **THEN** `hf_hub_download` SHALL return the cached path without re-downloading the full file

### Requirement: Cache directory excluded from git
The `.gitignore` at the repository root SHALL contain an entry that excludes `compare/VAD/.model_cache/` from version control.

#### Scenario: Model binaries not tracked
- **WHEN** a model is cached in `compare/VAD/.model_cache/`
- **THEN** `git status` SHALL NOT show any files from that directory as untracked or modified
