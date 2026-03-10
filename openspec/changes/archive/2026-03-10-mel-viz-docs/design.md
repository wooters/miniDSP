## Context

miniDSP's documentation is Doxygen-generated from markdown guide pages in `guides/`, with `guides/tutorials.md` serving as the index. The README mirrors the structure with "What's in the box?" feature lists and "Quick examples" code snippets. Every example in `examples/` has a corresponding guide page and README mention.

The mel_viz tool lives in `tools/mel_viz/` — a new directory convention distinct from `examples/`. It's not a simple API demo but a substantial program with its own `web/` sub-directory, Makefile, and README.

## Goals / Non-Goals

**Goals:**
- Add mel_viz to the documentation so users discover it
- Establish a "Tools" section pattern for future tools
- Keep the guide page lightweight — mel_viz already has its own `tools/mel_viz/README.md` with full usage docs

**Non-Goals:**
- Exhaustive tutorial-style walkthrough (mel_viz is a tool, not a DSP concept to teach)
- Adding mel_viz web assets to the Doxygen `HTML_EXTRA_FILES` (they're browser-served, not doc assets)
- Restructuring existing tutorials or examples

## Decisions

### 1. Separate "Tools" section in tutorials.md

**Decision**: Add a new "Tools" heading after the tutorial list in `guides/tutorials.md`, with its own `\subpage` entry. Don't intermix tools with tutorials.

**Rationale**: Tools are conceptually different from tutorials — tutorials teach DSP concepts via small API examples, while tools are standalone programs built on the library. A separate section makes this distinction clear and sets the pattern for future tools.

### 2. Guide page scope: overview + usage, not tutorial

**Decision**: The `guides/mel-viz.md` page should be a concise overview (what it does, how to build, how to run, architecture diagram) rather than a step-by-step tutorial. Link to `tools/mel_viz/README.md` for detailed usage and CLI flags.

**Rationale**: mel_viz's complexity lives in the JS renderer, not in DSP concepts. A tutorial-style "Reading the formula in C" section doesn't fit here. The tool's own README already has comprehensive usage docs.

### 3. README "Tools" section placement

**Decision**: Add a "Tools" section in the README after the "Quick examples" section and before "Python Bindings". It should have a brief description, a usage snippet, and a screenshot or link to the live demo.

**Rationale**: Tools are a higher-level showcase than examples — they demonstrate what you can build with the library. Placing them after "Quick examples" gives a natural progression from API usage → full applications.

### 4. No Doxyfile INPUT changes needed

**Decision**: The `guides/` directory is already in Doxygen's `INPUT` path, so adding `guides/mel-viz.md` is sufficient. No need to add `tools/` to INPUT — the tool's source code doesn't need to appear in the API reference.

**Rationale**: mel_viz is a consumer of the library, not part of the library API. Its internal code doesn't need Doxygen extraction.

## Risks / Trade-offs

- **Guide page could go stale**: mel_viz is in active development; the guide should reference the tool's own README for details rather than duplicating content. Mitigation: keep the guide page high-level.
- **Screenshot in docs**: A static screenshot of the visualizer doesn't capture the animated experience. Mitigation: mention that it's an interactive browser-based visualization and encourage running it.
