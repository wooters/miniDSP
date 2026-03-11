## Context

The mel_viz tool (`tools/mel_viz/web/`) is a single-page browser app with a vertical stacked layout: canvas → audio player → mode bar → export bar → controls → debug panel. All content is centered in a single column (max-width 600px). On shorter viewports, the controls are pushed below the fold, requiring users to scroll away from the visualization to adjust settings.

The codebase is vanilla HTML/CSS/JS with no build tooling or framework. Controls are generated dynamically by `controls.js` into a `<div id="controls">` container. The layout is driven entirely by `style.css` with a flexbox column on `body`.

## Goals / Non-Goals

**Goals:**
- Move all controls (palette, rings, sliders), mode bar, export bar, and debug panel into a collapsible left side panel
- Allow toggling the panel open/closed with a single button
- Persist panel state across page reloads via `localStorage`
- Keep the visualization and audio player centered in the remaining viewport space
- Maintain the existing dark theme aesthetic

**Non-Goals:**
- Drag-to-resize panel width (fixed width is sufficient)
- Responsive/mobile-first layout (the tool targets desktop browsers)
- Changes to `controls.js` API or `renderer.js` — only the DOM structure and CSS change
- Right-side or bottom-docked panel alternatives

## Decisions

### 1. CSS-driven panel with `transform: translateX`

The panel will be a fixed-width `<aside>` (260px) positioned on the left. Collapse/expand is handled by toggling a CSS class that applies `transform: translateX(-260px)` to slide it off-screen. A CSS `transition` on `transform` provides smooth animation.

**Why not `display: none`?** Transform preserves DOM state (slider positions, select values) and enables animation. `display: none` would require re-rendering controls on expand and can't be animated.

**Why not `position: fixed`?** Fixed positioning would overlay the panel on top of the canvas. Instead, the panel participates in normal flow (or uses a CSS grid/flex layout) so the main content area shifts when the panel opens.

### 2. Flexbox row layout on body

Change `body` from `flex-direction: column` to `flex-direction: row`. The `<aside>` is the first flex child; a new `<main>` wrapper holds the canvas, audio player, and title. The `<main>` element uses `flex: 1` to fill remaining width and keeps the existing column-centered layout internally.

When the panel collapses, its width transitions to 0 (via `margin-left: -260px` or `transform` + `width: 0` on a wrapper) and `<main>` expands to fill the viewport.

### 3. Toggle button position

A small toggle button (chevron icon `‹`/`›` or hamburger) will be placed at the top-right edge of the side panel (or just outside it in the main area). It stays visible in both states so users can always find it. Pure CSS/HTML — no icon library.

### 4. localStorage for state persistence

On toggle, write `melVizPanelOpen: "true"/"false"` to `localStorage`. On page load, read it and apply the appropriate class to the panel before first paint (inline `<script>` in `<head>` or early in the module) to avoid a flash of wrong layout.

## Risks / Trade-offs

- **Panel width vs. canvas centering** — With a 260px panel open, the canvas centering shifts right. This is acceptable since the canvas has a fixed rendered size (600px max display) and the main area will still be wide enough on any reasonable desktop viewport. → Mitigation: set a sensible `min-width` on `<main>` and test at 1024px viewport width.
- **Existing `controls.js` DOM generation** — `initControls(container, callback)` appends rows to whatever container it receives. As long as we pass the side panel's inner container as `controlsEl`, no changes to `controls.js` are needed.
- **Export bar width** — The export bar (`#export-bar`) currently has `max-width: 600px` and a flex row layout with a progress bar and percentage text, styled for the bottom-stacked layout. Inside a 260px panel, the `max-width` should be removed and the progress bar layout may need to stack vertically or shrink. The progress bar uses `flex: 1` which will adapt, but the overall layout may feel cramped. → Mitigation: remove `max-width: 600px` from `#export-bar` and let it fill the panel width; verify the progress bar and button fit at 260px.
- **Debug panel width** — The debug panel has `max-width: 600px` and a fixed-width meter canvas (580px). Inside a 260px panel, it will need narrower styling. → Mitigation: make the debug meter canvas `width: 100%` and use CSS to constrain it within the panel.
