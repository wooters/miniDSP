## ADDED Requirements

### Requirement: Side panel contains all controls and mode UI
The page SHALL render a side panel on the left edge containing, in order: mode bar (File/Mic buttons + status text), export bar, control widgets (palette, rings, smoothing, bass, wobble, glow), and debug panel. The main content area SHALL contain only the title, canvas visualization, and audio player.

#### Scenario: Page loads with precomputed data
- **WHEN** the page loads with `MEL_VIZ_DATA` available
- **THEN** the side panel contains the File button, Mic button, status text, export bar, and all control widgets
- **AND** the main content area contains only the title (`h1`), canvas, and audio player

#### Scenario: Page loads in mic-only mode
- **WHEN** the page loads without `MEL_VIZ_DATA`
- **THEN** the side panel contains the Mic button, status text, and all control widgets
- **AND** the export bar and File button are hidden

### Requirement: Panel is collapsible via toggle button
The side panel SHALL have a toggle button that collapses and expands it. The toggle button SHALL remain visible in both collapsed and expanded states.

#### Scenario: User collapses the panel
- **WHEN** the user clicks the toggle button while the panel is open
- **THEN** the panel slides out of view to the left
- **AND** the main content area expands to fill the full viewport width
- **AND** the toggle button remains visible at the left edge

#### Scenario: User expands the panel
- **WHEN** the user clicks the toggle button while the panel is collapsed
- **THEN** the panel slides into view from the left
- **AND** the main content area shifts right to accommodate the panel width

### Requirement: Panel toggle animates smoothly
The panel expand/collapse transition SHALL be animated using CSS transitions with a duration between 200ms and 300ms.

#### Scenario: Collapse animation
- **WHEN** the panel transitions from open to closed
- **THEN** the panel slides left with a smooth CSS transition (not instant)
- **AND** the main content area width change is also animated

### Requirement: Panel state persists across page reloads
The collapsed/expanded state of the panel SHALL be persisted in `localStorage` and restored on page load.

#### Scenario: User closes panel and reloads page
- **WHEN** the user collapses the panel and reloads the page
- **THEN** the panel loads in the collapsed state

#### Scenario: User opens panel and reloads page
- **WHEN** the user expands the panel and reloads the page
- **THEN** the panel loads in the expanded state

#### Scenario: First visit with no stored preference
- **WHEN** the page loads with no `localStorage` entry for panel state
- **THEN** the panel defaults to the expanded (open) state

### Requirement: Main content area centers visualization
The canvas and audio player in the main content area SHALL be horizontally centered within the available space, regardless of whether the side panel is open or closed.

#### Scenario: Panel open
- **WHEN** the side panel is expanded
- **THEN** the canvas and audio player are centered within the remaining viewport width (viewport minus panel width)

#### Scenario: Panel closed
- **WHEN** the side panel is collapsed
- **THEN** the canvas and audio player are centered within the full viewport width

### Requirement: Panel has fixed width
The side panel SHALL have a fixed width (not resizable by the user). The panel width SHALL be sufficient to display all controls without horizontal scrolling.

#### Scenario: Panel displays controls without overflow
- **WHEN** the side panel is expanded
- **THEN** all controls (dropdowns, sliders, labels, value displays) are fully visible without horizontal scrolling
