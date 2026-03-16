## 1. Doxyfile Configuration

- [x] 1.1 Add `doxygen-awesome-css/doxygen-awesome-interactive-toc.js` to `HTML_EXTRA_FILES` in `Doxyfile`

## 2. Custom Header

- [x] 2.1 Add `<script>` tag for `doxygen-awesome-interactive-toc.js` in `doxygen-custom/header.html`
- [x] 2.2 Add `DoxygenAwesomeInteractiveToc.init()` call in the initialization script block

## 3. Verification

- [x] 3.1 Run `make docs` and confirm `docs/html/doxygen-awesome-interactive-toc.js` exists
- [x] 3.2 Open a tutorial page in a browser and confirm the right-side TOC panel renders
