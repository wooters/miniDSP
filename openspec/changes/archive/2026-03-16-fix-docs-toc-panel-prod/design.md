## Context

The docs site at wooters.github.io/miniDSP is built by a GitHub Actions workflow (`docs.yml`) on ubuntu-24.04. The `apt` package for Doxygen on that platform is version **1.9.8**. Locally, the project uses Doxygen **1.16.1** (Homebrew).

Doxygen 1.16 introduced a native right-side "On This Page" panel (`<div id="page-nav" class="page-nav-panel">`) populated dynamically by `navtree.js`. This panel does not exist in 1.9.8.

The custom `header.html` is templated for 1.16.1 (first line: `<!-- HTML header for doxygen 1.16.1-->`). Using it with 1.9.8 may cause additional template-variable mismatches beyond the missing panel.

The prior fix (`fix-docs-toc-panel`) wired up `doxygen-awesome-interactive-toc.js`, which is a separate mechanism requiring `[TOC]` markers in markdown. No guide pages have `[TOC]`, so the extension never activates — it's dead code.

## Goals / Non-Goals

**Goals:**
- Deployed docs at wooters.github.io/miniDSP show the right-side TOC panel, matching local builds.
- CI and local Doxygen versions are aligned (both 1.16.1).

**Non-Goals:**
- Adding `[TOC]` markers to guide pages (the native panel supersedes this approach).
- Upgrading to a bleeding-edge Doxygen version beyond 1.16.1.
- Changing Doxygen configuration settings (GENERATE_TREEVIEW, FULL_SIDEBAR, etc.).

## Decisions

### 1. Install Doxygen 1.16.1 from upstream binary on CI

**Rationale:** The ubuntu-24.04 apt package (1.9.8) is too old. Doxygen publishes official Linux binaries on GitHub Releases. Download and extract the binary in the workflow instead of using `apt install doxygen`.

**Alternatives considered:**
- *Use a PPA*: No official Doxygen PPA with 1.16.x. Third-party PPAs are a maintenance risk.
- *Build from source on CI*: Slow (~5 min), adds cmake/flex/bison dependencies. Unnecessary when binaries exist.
- *Pin a Docker image with Doxygen 1.16*: Adds container complexity for one tool.

### 2. Remove dead doxygen-awesome-interactive-toc wiring

**Rationale:** The extension requires `.contents > .toc` DOM elements, generated only when markdown has `[TOC]` markers. No guide pages use `[TOC]`, so the JS loads, calls `init()`, finds nothing, and exits. The native 1.16 panel is the actual right-side TOC. Keeping dead wiring is confusing for future maintainers.

**What to remove:**
- `doxygen-custom/header.html`: Remove the `<script>` tag for `doxygen-awesome-interactive-toc.js` and the `DoxygenAwesomeInteractiveToc.init()` call.
- `Doxyfile`: Remove `doxygen-awesome-interactive-toc.js` from `HTML_EXTRA_FILES`.

## Risks / Trade-offs

- **Binary URL stability** → Doxygen GitHub Releases have stable URLs. Pin the exact version (1.16.1) so breakage is caught immediately if the URL changes.
- **CI build time** → Downloading a ~5 MB tarball adds <5 seconds. Negligible.
- **Future Doxygen upgrades** → Pinning the version means manual bumps. This is intentional — Doxygen version changes can alter HTML output and break custom headers.
- **header.html template drift** → If Doxygen 1.16.1 is later replaced, the custom header must be regenerated (`doxygen -w html`). Document this in CLAUDE.md.
