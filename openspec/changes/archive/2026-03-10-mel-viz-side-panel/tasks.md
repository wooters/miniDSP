## 1. Restructure HTML

- [x] 1.1 Wrap canvas, title, and audio player in a `<main>` element
- [x] 1.2 Create an `<aside id="side-panel">` containing: toggle button, mode bar, export bar, controls div, and debug panel
- [x] 1.3 Add a toggle button element (chevron `‹`/`›`) inside the aside at the top-right edge

## 2. Rework CSS layout

- [x] 2.1 Change `body` from `flex-direction: column` to `flex-direction: row`
- [x] 2.2 Style the `<aside>` with fixed width (260px), dark background, overflow-y auto, full viewport height
- [x] 2.3 Style `<main>` with `flex: 1`, column layout, centered content
- [x] 2.4 Add `.collapsed` class on `<aside>` that uses `margin-left: -260px` (or transform) to slide it off-screen
- [x] 2.5 Add CSS transition (250ms) on the aside margin/transform for smooth animation
- [x] 2.6 Style the toggle button to remain visible in both states (positioned at the panel edge)
- [x] 2.7 Adapt export bar styles (remove `max-width: 600px`, verify progress bar + button fit) for 260px panel
- [x] 2.8 Adapt debug panel styles (meter canvas width, max-width) to fit within 260px panel

## 3. Wire up toggle behavior

- [x] 3.1 Add click handler on the toggle button to add/remove `.collapsed` class on the aside
- [x] 3.2 Update toggle button icon (chevron direction) on state change
- [x] 3.3 Read `localStorage("melVizPanelOpen")` on page load and apply initial collapsed/expanded state
- [x] 3.4 Write panel state to `localStorage` on each toggle

## 4. Verify and clean up

- [x] 4.1 Confirm `initControls(controlsEl, ...)` still works with the controls div inside the aside
- [x] 4.2 Test file mode and mic mode with panel open and closed
- [x] 4.3 Remove any leftover bottom-layout CSS rules that no longer apply
