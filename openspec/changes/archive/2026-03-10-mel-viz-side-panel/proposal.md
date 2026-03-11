## Why

The mel visualizer's controls (palette, rings, smoothing, bass, wobble, glow) currently sit below the canvas and audio player, pushing them off-screen on shorter viewports and requiring scrolling to adjust settings while watching the visualization. Moving the controls to a collapsible side panel on the left keeps them accessible without competing for vertical space, and lets users hide them entirely to maximize the visualization area.

## What Changes

- Replace the current bottom-stacked controls layout with a collapsible side panel on the left edge of the page.
- Add a toggle button (hamburger/chevron) to expand/collapse the panel.
- When the panel is collapsed, the visualization and audio player center in the full viewport width.
- When expanded, the panel slides in from the left and the main content area shifts right to accommodate it.
- Move the mode bar (File/Mic buttons + status), export bar, and debug panel into the side panel alongside the existing controls.
- Persist collapsed/expanded state in `localStorage` so the preference survives page reloads.

## Capabilities

### New Capabilities
- `side-panel-layout`: Collapsible left side panel containing all controls, mode switching, and export UI. Covers panel toggle behavior, layout responsiveness, and state persistence.

### Modified Capabilities

(none — no existing spec-level requirements are changing)

## Impact

- **Files changed**: `index.html` (DOM restructure), `style.css` (layout overhaul from vertical to horizontal with panel), `controls.js` (minor — may need to render into a new container structure)
- **No API or library changes** — this is purely a front-end layout change within `tools/mel_viz/web/`.
- **No new dependencies** — CSS transitions/transforms handle the collapse animation.
