## ADDED Requirements

### Requirement: CI installs Doxygen 1.16.1 from upstream binary
The docs workflow SHALL install Doxygen 1.16.1 from the official Doxygen GitHub Releases binary (Linux x86_64) instead of the ubuntu-24.04 apt package.

#### Scenario: Docs workflow builds with Doxygen 1.16.1
- **WHEN** the docs workflow runs on CI
- **THEN** `doxygen --version` reports `1.16.1`

#### Scenario: Generated HTML contains page-nav panel
- **WHEN** `make docs` completes on CI
- **THEN** guide pages (e.g., `signal-generators.html`) contain `<div id="page-nav" class="page-nav-panel">`

### Requirement: Dead interactive-toc wiring is removed
The `doxygen-awesome-interactive-toc.js` script tag, `init()` call, and `HTML_EXTRA_FILES` entry SHALL be removed since no pages use `[TOC]` markers and the native 1.16 panel supersedes this extension.

#### Scenario: header.html does not load interactive-toc
- **WHEN** inspecting `doxygen-custom/header.html`
- **THEN** there is no `<script>` tag referencing `doxygen-awesome-interactive-toc.js` and no `DoxygenAwesomeInteractiveToc.init()` call

#### Scenario: Doxyfile does not list interactive-toc in EXTRA_FILES
- **WHEN** inspecting `Doxyfile`
- **THEN** `HTML_EXTRA_FILES` does not contain `doxygen-awesome-interactive-toc.js`
