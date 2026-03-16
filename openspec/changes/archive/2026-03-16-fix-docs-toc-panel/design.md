## Context

The Doxygen Awesome theme (v2.4.1 submodule) ships with an interactive TOC extension (`doxygen-awesome-interactive-toc.js`) that renders a right-side page outline on content pages. The CSS for this feature is already included via `doxygen-awesome.css`, but the JS file is neither copied to the output directory nor loaded in the custom header. The other three extensions (dark mode toggle, fragment copy button, paragraph link) are correctly wired up — the TOC was simply missed.

## Goals / Non-Goals

**Goals:**
- Enable the interactive TOC panel on all docs pages by wiring up the existing JS extension.

**Non-Goals:**
- Customizing the TOC behavior or appearance beyond what the extension provides out of the box.
- Changing the Doxygen Awesome submodule version.

## Decisions

**Wire up the existing extension following the same pattern as the other three.**

The three active extensions follow a consistent pattern:
1. Listed in `HTML_EXTRA_FILES` in the Doxyfile
2. Loaded via `<script>` tag in `doxygen-custom/header.html`
3. Initialized via `ExtensionName.init()` call in the same script block

The TOC extension will follow the same pattern. No alternative approaches are needed — this is purely a configuration gap.

## Risks / Trade-offs

- **[Low] Layout shift on narrow viewports** — The TOC collapses into a mobile-friendly menu on screens < 1000px (handled by the extension's built-in CSS). No custom responsive work needed.
- **[Low] Pages with no headings** — The extension gracefully hides when there are no TOC-eligible headings. No edge-case handling needed.
