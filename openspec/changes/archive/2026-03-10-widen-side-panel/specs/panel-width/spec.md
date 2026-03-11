## ADDED Requirements

### Requirement: Side panel width provides readable slider controls
The side panel SHALL be 340px wide, providing sufficient space for slider tracks to be easily adjusted and their values read.

#### Scenario: Panel displays at correct width
- **WHEN** the mel_viz page loads with the side panel expanded
- **THEN** the side panel SHALL be 340px wide

#### Scenario: Panel collapse still works
- **WHEN** the user clicks the collapse toggle
- **THEN** the panel SHALL slide fully off-screen using a negative margin equal to the panel width
