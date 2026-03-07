## Why

The `pyminidsp` Python wrapper (wooters/pyminidsp) needs a way to pin to known-good versions of the C library. Without git tags or version identifiers, the wrapper must reference raw commit SHAs. Adding semver tags and a `make install` target gives downstream consumers a stable, conventional integration point.

## What Changes

- Add a `VERSION` file at the repo root containing the current semver string (starting at `0.1.0`)
- Add a `MINIDSP_VERSION` macro to `include/minidsp.h` that embeds the version at compile time
- Add a `make install` target that copies headers and the static library to a configurable prefix
- Add a `make uninstall` target for clean removal
- Create an initial `v0.1.0` git tag

## Capabilities

### New Capabilities
- `library-versioning`: Version file, header macro, and git tag conventions
- `make-install`: Install/uninstall targets for the static library and headers

### Modified Capabilities
<!-- None -->

## Impact

- **Build system**: New `install`/`uninstall` targets in root Makefile; new `VERSION` file
- **Public headers**: `minidsp.h` gains `MINIDSP_VERSION` / `MINIDSP_VERSION_MAJOR` / `MINIDSP_VERSION_MINOR` / `MINIDSP_VERSION_PATCH` macros
- **Git workflow**: Releases tagged with `vMAJOR.MINOR.PATCH`
- **Downstream (pyminidsp)**: Can pin via git tag, read version from header or file, and use `make install PREFIX=...` for system-level builds
