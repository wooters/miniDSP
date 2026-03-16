## Why

The right-side "On This Page" TOC panel works locally (Doxygen 1.16.1) but is missing on the deployed GitHub Pages site. The CI workflow (`docs.yml`) installs Doxygen from ubuntu-24.04 apt, which ships version 1.9.8. The native `page-nav-panel` feature was introduced in Doxygen 1.16+ — it simply doesn't exist in 1.9.8.

## What Changes

- Upgrade the CI docs workflow to install Doxygen 1.16.1 (matching local) instead of the ubuntu-24.04 apt package (1.9.8).
- The custom `header.html` template is already written for 1.16.1 — no changes needed there.
- Remove the `doxygen-awesome-interactive-toc.js` wiring (script tag + `init()` call in header.html and `HTML_EXTRA_FILES` entry in Doxyfile) since it's dead code — no guide pages use `[TOC]` markers, so the extension never activates. The native 1.16.1 panel supersedes it.

## Capabilities

### New Capabilities
- `ci-doxygen-upgrade`: Upgrade the GitHub Actions docs workflow to install Doxygen 1.16.1 from upstream binary instead of the outdated apt package.

### Modified Capabilities

(none)

## Impact

- `.github/workflows/docs.yml` — changed install step
- `doxygen-custom/header.html` — remove dead interactive-toc script/init lines
- `Doxyfile` — remove `doxygen-awesome-interactive-toc.js` from `HTML_EXTRA_FILES`
- No API, library code, or test changes
