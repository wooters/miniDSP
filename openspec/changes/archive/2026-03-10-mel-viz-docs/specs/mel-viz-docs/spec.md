# mel-viz-docs

Documentation for the mel_viz audio visualizer tool.

## Requirements

### Guide page (`guides/mel-viz.md`)
- Page anchor: `{#mel-viz}`
- Sections: what it does, build, file mode usage, mic mode usage, architecture overview
- Include a simple ASCII or text architecture diagram (C program → data.js → browser renderer)
- Link to `tools/mel_viz/README.md` for detailed CLI flags and knob descriptions
- Do NOT include tutorial-style "Reading the formula in C" sections (this is a tool overview, not a DSP tutorial)

### Tutorials index (`guides/tutorials.md`)
- Add a `## Tools` heading after the existing tutorial list
- Add `\subpage mel-viz` entry under the Tools heading
- Brief description: mel-spectrum audio visualizer with browser-based radial animation

### README (`README.md`)
- Add a `## Tools` section after "Quick examples" and before "Python Bindings"
- Include mel_viz with: one-sentence description, build command, usage snippet, link to guide page
- Keep it concise — 10-15 lines max

### Doxyfile
- No changes expected (guides/ is already in INPUT)
- If the guide references snippets from `tools/`, add `tools/mel_viz` to `EXAMPLE_PATH`
