## Context

The README is the project's front door. Four issues have been identified: no platform info, a missing tool listing, suboptimal section order, and an incorrect dependency note.

## Goals / Non-Goals

**Goals:**
- Make the README accurate and complete for someone encountering the project for the first time
- Ensure all tools in `tools/` are documented

**Non-Goals:**
- Rewriting prose or restructuring beyond the four requested changes
- Adding Windows build support (just noting it's untested)

## Decisions

1. **Platform note placement** — Add a short sentence right after the project description (line 8 area), before the "What's in the box?" section. This is the natural place a reader looks for compatibility info.

2. **audio_steg tool description** — Mirror the style of the existing `mel_viz` entry: a `###` heading, one-line description, build/run snippet, and link to docs (if available). Source `tools/audio_steg/audio_steg.c` header comment for the description.

3. **Section reorder** — Move the entire "Build and Test" section (with all subsections: Dependencies, Compile, Test, Container, Docs, Hooks) to appear before "Use in your project". Rationale: readers need to install deps and build before they can use the library.

4. **Apple container fix** — Change "macOS 26+ built-in" to something like "Install from [GitHub](https://github.com/apple/container)" to accurately reflect that it's a separate download.

## Risks / Trade-offs

- **Link rot** — The audio_steg docs link may not exist yet on GitHub Pages. Use a conditional link or omit until docs are published.
- **Section reorder affects anchor links** — No external links to README anchors are known, so this is low risk.
