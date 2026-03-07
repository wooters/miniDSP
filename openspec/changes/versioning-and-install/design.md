## Context

miniDSP is a static C library with no versioning scheme. The new `pyminidsp` Python wrapper needs to reference stable snapshots. The project uses a shared `config.mk`, a root Makefile, and Homebrew-based dependency resolution on macOS.

## Goals / Non-Goals

**Goals:**
- Single source of truth for the version number
- Compile-time version macros for downstream C consumers
- Standard `make install`/`uninstall` for system-level integration
- Git tags following semver convention

**Non-Goals:**
- Shared library (`.so`/`.dylib`) support — static library only for now
- Automated release workflows or CI-based tagging
- pkg-config `.pc` file generation (can be added later)

## Decisions

### 1. Single source of truth: `VERSION` file

**Decision**: A plain-text `VERSION` file at the repo root containing just the semver string (e.g., `0.1.0`).

**Rationale**: Simple to parse from shell (`cat VERSION`), Python (`open("VERSION").read().strip()`), and Makefile (`$(shell cat VERSION)`). No build step needed to extract the version.

**Alternative considered**: Encoding version only in the header `#define` — harder to parse from non-C tooling (pyminidsp's build system, release scripts).

### 2. Header macros derived from `VERSION` file

**Decision**: The Makefile reads `VERSION`, splits it into major/minor/patch, and passes them as `-D` flags. `minidsp.h` provides:
- `MINIDSP_VERSION_MAJOR`, `MINIDSP_VERSION_MINOR`, `MINIDSP_VERSION_PATCH` (integers)
- `MINIDSP_VERSION` (string literal, e.g., `"0.1.0"`)

These are defined via `-D` in CFLAGS so the `VERSION` file remains the single source of truth. The header provides `#ifndef` guards with defaults so the header is still self-contained if someone compiles without the Makefile.

**Rationale**: Downstream C code can do compile-time version checks. The `-D` approach avoids needing a generated header file.

### 3. Install target layout

**Decision**: `make install` copies to a configurable `PREFIX` (default `/usr/local`):
```
$(PREFIX)/lib/libminidsp.a
$(PREFIX)/include/minidsp.h
$(PREFIX)/include/biquad.h
$(PREFIX)/include/fileio.h
$(PREFIX)/include/liveio.h
```

**Rationale**: Follows the standard Unix convention. `PREFIX` is overridable (`make install PREFIX=/opt/minidsp`). Only public headers from `include/` are installed — internal headers from `src/` are not.

### 4. Starting version: `0.1.0`

**Decision**: Start at `0.1.0` (pre-1.0 semver) to signal that the API is not yet stable.

**Rationale**: The library is functional but the API may still evolve. `0.x.y` communicates this clearly per semver conventions.

### 5. Tag format: `vMAJOR.MINOR.PATCH`

**Decision**: Annotated git tags with `v` prefix (e.g., `v0.1.0`).

**Rationale**: The `v` prefix is the most common convention for C projects and is what GitHub recognizes for release pages.

## Risks / Trade-offs

- **[Low] VERSION file gets out of sync with git tag** → Mitigated by documenting the release checklist: update VERSION, commit, tag, push.
- **[Low] Install target conflicts with system packages** → Default prefix `/usr/local` is standard for user-installed software. Users can override with `PREFIX=`.
