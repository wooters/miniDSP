## ADDED Requirements

### Requirement: Brickwall lowpass guide page exists

A guide page SHALL exist at `guides/brickwall-lowpass.md` with Doxygen page anchor `{#brickwall-lowpass}`. The guide SHALL contain:

1. An introductory paragraph explaining FFT-based brickwall lowpass filtering with a Wikipedia link to [lowpass filters](https://en.wikipedia.org/wiki/Low-pass_filter)
2. The frequency-domain filter formula using `\f[...\f]` display math (matching the formula in `minidsp.h`)
3. A "Reading the formula in C" section showing the bin-zeroing loop with explicit variable mapping comments
4. An API section with a quick code example calling `MD_lowpass_brickwall()`
5. A note on Gibbs ringing as the primary tradeoff
6. An embedded example output (iframe) showing before/after magnitude spectra

#### Scenario: Guide page renders in Doxygen
- **WHEN** `make docs` is run
- **THEN** the brickwall lowpass guide page SHALL appear in the generated HTML documentation with rendered formulas, code blocks, and embedded iframe

### Requirement: Tutorials index includes brickwall guide

`guides/tutorials.md` SHALL contain a `\subpage brickwall-lowpass` entry. The entry SHALL be positioned after the FIR convolution entry and before the pitch detection entry.

#### Scenario: Tutorials index ordering
- **WHEN** a user views the tutorials index page
- **THEN** the brickwall lowpass entry SHALL appear between "FIR Filters and Convolution" and "Pitch Detection"

### Requirement: Example program exists

An example program SHALL exist at `examples/brickwall.c` that:
1. Generates a mixed signal with a low-frequency and high-frequency tone
2. Computes the magnitude spectrum before filtering
3. Applies `MD_lowpass_brickwall()` at a chosen cutoff
4. Computes the magnitude spectrum after filtering
5. Writes a self-contained HTML file with a two-subplot Plotly chart showing before/after spectra
6. Contains `//! [snippet-name]` markers for embedding in the guide via `\snippet`

#### Scenario: Example compiles and runs
- **WHEN** `make -C examples brickwall` is run
- **AND** `./examples/brickwall` is executed
- **THEN** it SHALL produce `examples/brickwall.html` with a valid Plotly chart

### Requirement: Example integrated into build

The example program SHALL be added to:
1. `examples/Makefile` `EXAMPLES` variable
2. `examples/Makefile` `plot` target
3. `.gitignore` with lines for `examples/brickwall`, `examples/brickwall.csv`, `examples/brickwall.html`
4. `.dockerignore` with the same three lines

#### Scenario: Clean build includes brickwall
- **WHEN** `make -C examples` is run from a clean state
- **THEN** the `brickwall` binary SHALL be built alongside other examples

### Requirement: Doxyfile updated for iframe asset

`Doxyfile` `HTML_EXTRA_FILES` SHALL include the brickwall example HTML output so it can be embedded as an iframe in the guide.

#### Scenario: Iframe source available in docs
- **WHEN** `make docs` completes (including running the example to generate the HTML)
- **THEN** the brickwall HTML file SHALL be copied to `docs/html/` and accessible as an iframe source
