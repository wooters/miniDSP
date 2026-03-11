## Why

The mel_viz side panel is fixed at 260px wide, which leaves only ~96px for slider tracks after accounting for labels (80px), value displays (40px), gaps, and padding. This makes it difficult to read slider values and precisely adjust controls like Smoothing, Bass, Wobble, and Glow.

## What Changes

- Increase the side panel width from 260px to 340px
- Update the collapse margin to match the new width
- Ensure slider tracks have adequate space for fine-grained control

## Capabilities

### New Capabilities

### Modified Capabilities

## Impact

- `tools/mel_viz/web/style.css` — panel width and collapse margin values
