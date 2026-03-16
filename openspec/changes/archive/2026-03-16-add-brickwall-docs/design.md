## Context

`MD_lowpass_brickwall()` is implemented and tested but has no tutorial guide page. The existing guides follow a consistent pattern: theory with LaTeX formulas, "Reading the formula in C" walkthroughs, API reference, a working example program with snippet markers, and embedded HTML visualizations. The function is currently documented only in its header doc-comment and referenced in the audio steganography pipeline diagram.

## Goals / Non-Goals

**Goals:**
- Add a standalone tutorial guide page that teaches the brickwall lowpass concept and API
- Follow the established guide conventions (formulas, C walkthroughs, snippets, plots)
- Place the guide logically near existing filtering content in the tutorials index

**Non-Goals:**
- Covering general IIR/biquad filter design (separate topic)
- Adding interactive or animated visualizations (a static before/after magnitude spectrum is sufficient)
- Modifying the library implementation or API

## Decisions

### 1. Guide placement in tutorials index

**Decision**: Insert the `\subpage` entry after the FIR convolution entry and before the pitch detection entry.

- *Alternative: Group with resampling*: Resampling is about sample rate conversion, not filtering per se. The brickwall is a frequency-domain filter, most naturally grouped with FIR convolution.
- *Alternative: Append at end*: Would break the thematic flow — filtering topics should be adjacent.

**Rationale**: The FIR convolution guide already covers time-domain lowpass filters. The brickwall guide extends this with the frequency-domain approach, forming a natural progression.

### 2. Example program design

**Decision**: Create `examples/brickwall.c` that generates a mixed signal (low + high frequency tones), applies `MD_lowpass_brickwall()`, and writes a before/after magnitude spectrum comparison to an HTML file using Plotly (two-subplot chart).

- *Alternative: CSV output with gnuplot*: Most recent examples use self-contained HTML with Plotly, matching the project's newer convention.
- *Alternative: Time-domain plot*: The frequency domain is where the brickwall effect is most visible and instructive.

**Rationale**: A before/after spectrum plot is the clearest way to show perfect suppression above the cutoff. Plotly HTML is self-contained and consistent with recent examples.

### 3. Guide structure

**Decision**: Follow the single-page pattern used by other spectrum-related guides (magnitude-spectrum, power-spectral-density):

1. Intro paragraph with Wikipedia link
2. Formula with `\f[...\f]` display math
3. "Reading the formula in C" section
4. API section with quick example
5. Gibbs ringing note (tradeoff awareness)
6. Embedded example output (iframe)

**Rationale**: Matches established guide conventions exactly. The Gibbs ringing note is important because it's the main practical tradeoff of brickwall filters.

## Risks / Trade-offs

**[Guide may feel thin]** → The brickwall filter is conceptually simple (zero bins above cutoff). This is actually a feature — the guide serves as a concise, self-contained reference. No need to pad it.

**[Iframe HTML asset management]** → The example HTML output needs to be added to `Doxyfile` `HTML_EXTRA_FILES` for iframe embedding, or the guide can use a `@code` block showing the plot inline. Decision: use `HTML_EXTRA_FILES` to match the existing pattern.
