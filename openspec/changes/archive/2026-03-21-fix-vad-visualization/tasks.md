## 1. Rename visualization output

- [x] 1.1 In `examples/vad.c`, change the output filename from `"vad.html"` to `"vad_plot.html"`
- [x] 1.2 In `guides/vad.md`, update the iframe `src` from `"vad.html"` to `"vad_plot.html"`
- [x] 1.3 In `Doxyfile`, update `HTML_EXTRA_FILES` entry from `examples/vad.html` to `examples/vad_plot.html`

## 2. Update ignore patterns

- [x] 2.1 In `.gitignore`, update `examples/vad.html` to `examples/vad_plot.html`
- [x] 2.2 In `.dockerignore`, update `examples/vad.html` to `examples/vad_plot.html`

## 3. Verify

- [x] 3.1 Build the example (`make -C examples vad`) and confirm `examples/vad_plot.html` is generated
- [x] 3.2 Build docs (`make docs`) and confirm the guide page and visualization are separate files in `docs/html/`
