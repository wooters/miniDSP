## 1. Guide page

- [x] 1.1 Create `guides/mel-viz.md` with `{#mel-viz}` anchor — overview, build, file mode, mic mode, architecture diagram
- [x] 1.2 Link to `tools/mel_viz/README.md` for detailed CLI flags and visual knob descriptions
- [x] 1.3 Mention `samples/` directory contains audio files (e.g., `punchy_slap_bass_30s.wav`) that can be used with mel_viz

## 2. Tutorials index

- [x] 2.1 Add `## Tools` heading at the end of `guides/tutorials.md`
- [x] 2.2 Add `\subpage mel-viz` entry with brief description under the Tools heading

## 3. README

- [x] 3.1 Add `## Tools` section in `README.md` after "Quick examples" and before "Python Bindings"
- [x] 3.2 Add mel_viz entry with description, build command (`make tools`), usage snippet, and link to docs page
- [x] 3.3 Mention `samples/` directory in the README Tools section or near the mel_viz usage snippet

## 4. Build verification

- [x] 4.1 Run `make docs` and verify the mel-viz guide page renders correctly in the generated HTML
- [x] 4.2 Verify the tutorials index links to the new page
