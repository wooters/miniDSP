## Context

The codebase was migrated from C23 to C17 (commit `9933eca`) for broader compiler compatibility, but documentation artifacts (README, Dockerfile, badge) were not fully updated. The build system (`config.mk`) correctly uses `-std=c17`, so this is purely a docs/config alignment fix.

Current state of stale references:
- `README.md` line 6: badge says "C23"
- `README.md` line 146, 156: compile examples use `-std=c23`
- `README.md` line 439: states GCC 14 required for `-std=c23`
- `README.md` line 455: container description mentions "GCC 14"
- `Dockerfile` line 4: explicitly installs `gcc-14` and sets up alternatives

## Goals / Non-Goals

**Goals:**
- All documentation accurately reflects the C17 standard and its compiler requirements
- Dockerfile uses the simplest GCC setup that works with C17
- Users on Ubuntu 24.04 (default GCC 13) can build without installing extra packages

**Non-Goals:**
- Changing the actual C standard used by the build system (already C17)
- Modifying any library code or tests
- Changing CI workflows beyond what's needed for the Dockerfile update

## Decisions

**Use system default `gcc` in Dockerfile instead of `gcc-14`:**
C17 is fully supported by GCC 8+. Ubuntu 24.04 ships GCC 13, which is more than sufficient. Removing the explicit `gcc-14` install simplifies the Dockerfile and eliminates the `update-alternatives` step. The `gcc` package on Ubuntu 24.04 provides GCC 13 via the default `gcc` command.

**Remove the GCC 14 requirement paragraph entirely:**
Rather than softening it to "GCC 13 or later", remove the paragraph since C17 support is ubiquitous in any modern GCC. No special compiler installation guidance is needed.

**Update badge to link to C17 Wikipedia page:**
The badge currently links to the C23 Wikipedia article. Update both the label and the URL.

## Risks / Trade-offs

**Risk: Dockerfile `gcc` package may differ across Ubuntu versions** → The Dockerfile is pinned to `ubuntu:24.04`, so the default `gcc` package is deterministic (GCC 13.2). If the base image is updated in the future, C17 will still be supported by any version they ship.

**Risk: Existing container images cached with gcc-14** → No mitigation needed; `make container-test` rebuilds the image. Users pulling fresh will get the simplified Dockerfile.
