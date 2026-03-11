## Context

The mel_viz side panel (`#side-panel` in `style.css`) is fixed at 260px. After padding (32px), label (80px), gaps (24px), and value display (40px), only ~96px remains for slider tracks. This makes fine adjustment difficult and values hard to read at a glance.

## Goals / Non-Goals

**Goals:**
- Widen the panel so slider tracks have ~170px+ of usable space
- Maintain collapse/expand toggle behavior

**Non-Goals:**
- Responsive/adaptive panel sizing
- Changing label widths or control layout
- Adding new controls

## Decisions

**Panel width: 340px** — Increases slider track space from ~96px to ~176px (nearly double). Still leaves the majority of the viewport for the visualization canvas. 340px is a common sidebar width in web apps.

**Approach: CSS-only change** — Update `width` and `margin-left` (collapse offset) in `style.css`. No JS changes needed since the collapse toggle already reads the panel width dynamically via the existing `margin-left` animation.

## Risks / Trade-offs

- [Slightly less canvas space] → 80px narrower visualization area. Acceptable since the canvas is typically >1000px wide on modern displays.
