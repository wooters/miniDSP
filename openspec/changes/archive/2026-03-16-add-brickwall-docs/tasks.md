## 1. Example Program

- [x] 1.1 Create `examples/brickwall.c`: generate a mixed signal (e.g., 400 Hz + 12000 Hz at 48 kHz), compute magnitude spectrum before and after `MD_lowpass_brickwall()`, write a two-subplot Plotly HTML file (`brickwall.html`). Include `//! [snippet-name]` markers for guide embedding.
- [x] 1.2 Add `brickwall` to `EXAMPLES` variable in `examples/Makefile`
- [x] 1.3 Add `./brickwall` run line to the `plot` target in `examples/Makefile`
- [x] 1.4 Add `examples/brickwall`, `examples/brickwall.csv`, `examples/brickwall.html` lines to `.gitignore`
- [x] 1.5 Add `examples/brickwall`, `examples/brickwall.csv`, `examples/brickwall.html` lines to `.dockerignore`

## 2. Guide Page

- [x] 2.1 Create `guides/brickwall-lowpass.md` with `{#brickwall-lowpass}` anchor: intro paragraph, display-math formula, "Reading the formula in C" section, API quick example, Gibbs ringing note, and embedded iframe showing example output
- [x] 2.2 Add `\subpage brickwall-lowpass` entry to `guides/tutorials.md` after the FIR convolution line and before the pitch detection line

## 3. Doxyfile

- [x] 3.1 Add the brickwall example HTML output path to `HTML_EXTRA_FILES` in `Doxyfile`

## 4. Verification

- [x] 4.1 Run `make -C examples brickwall` and execute `./examples/brickwall` to confirm HTML output is generated
- [x] 4.2 Run `make docs` and verify the brickwall guide page renders with formulas, code, and iframe
