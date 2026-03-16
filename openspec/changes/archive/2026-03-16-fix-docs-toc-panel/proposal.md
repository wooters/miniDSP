## Why

The Doxygen Awesome theme includes an interactive table-of-contents (TOC) panel that appears on the right side of content pages. The JS file (`doxygen-awesome-interactive-toc.js`) exists in the submodule but is not wired up in the Doxyfile or custom header, so the TOC panel never renders in the deployed docs on `wooters.github.io`.

## What Changes

- Add `doxygen-awesome-interactive-toc.js` to `HTML_EXTRA_FILES` in the Doxyfile so Doxygen copies it to the output directory.
- Add a `<script>` tag and `DoxygenAwesomeInteractiveToc.init()` call in the custom header so the TOC JS loads and activates on every page.

## Capabilities

### New Capabilities
- `interactive-toc`: Wire up the existing Doxygen Awesome interactive TOC extension so it renders a right-side page outline on tutorial and API pages.

### Modified Capabilities

(none)

## Impact

- **Files changed**: `Doxyfile`, `doxygen-custom/header.html`
- **No code changes** to the C library, tests, or build system.
- **No new dependencies** — the JS file already exists in the `doxygen-awesome-css` submodule.
