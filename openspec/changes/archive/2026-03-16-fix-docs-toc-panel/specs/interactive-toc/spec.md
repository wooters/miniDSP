## ADDED Requirements

### Requirement: Interactive TOC JS is copied to docs output
The Doxyfile SHALL list `doxygen-awesome-css/doxygen-awesome-interactive-toc.js` in `HTML_EXTRA_FILES` so Doxygen copies it to the output directory during generation.

#### Scenario: JS file present in output
- **WHEN** `make docs` is run
- **THEN** `docs/html/doxygen-awesome-interactive-toc.js` SHALL exist

### Requirement: Interactive TOC JS is loaded and initialized
The custom HTML header SHALL include a `<script>` tag loading `doxygen-awesome-interactive-toc.js` and SHALL call `DoxygenAwesomeInteractiveToc.init()` to activate the extension.

#### Scenario: TOC renders on a tutorial page with headings
- **WHEN** a tutorial page with multiple `##` headings is viewed in a browser
- **THEN** a right-side table of contents panel SHALL be visible listing those headings

#### Scenario: TOC hidden on pages without headings
- **WHEN** a page with no TOC-eligible headings is viewed
- **THEN** no TOC panel SHALL be displayed (extension hides itself gracefully)
