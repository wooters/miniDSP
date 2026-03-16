## 1. Upgrade Doxygen on CI

- [x] 1.1 Update `.github/workflows/docs.yml` to download and install Doxygen 1.16.1 Linux x86_64 binary from the official GitHub Releases, replacing the `apt install doxygen` approach
- [x] 1.2 Verify the workflow still checks out submodules and installs other build deps (gcc-14, libfftw3-dev, libsndfile1-dev, etc.)

## 2. Remove dead interactive-toc wiring

- [x] 2.1 Remove the `<script>` tag for `doxygen-awesome-interactive-toc.js` from `doxygen-custom/header.html`
- [x] 2.2 Remove the `DoxygenAwesomeInteractiveToc.init()` call from `doxygen-custom/header.html`
- [x] 2.3 Remove `doxygen-awesome-css/doxygen-awesome-interactive-toc.js` from `HTML_EXTRA_FILES` in `Doxyfile`

## 3. Verify

- [x] 3.1 Run `make docs` locally and confirm `docs/html/signal-generators.html` still contains `<div id="page-nav" class="page-nav-panel">`
- [ ] 3.2 Push to main and confirm the deployed site at wooters.github.io/miniDSP shows the right-side TOC panel on guide pages
