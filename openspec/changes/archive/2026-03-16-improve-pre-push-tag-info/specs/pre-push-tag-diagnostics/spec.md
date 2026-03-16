## ADDED Requirements

### Requirement: Display tag status header
The pre-push hook SHALL display a clearly delimited "VERSION / TAG STATUS" section after all tests pass when pushing to main.

#### Scenario: Tests pass and tag status section appears
- **WHEN** all tests pass on a push to main
- **THEN** the hook prints a section bounded by `========` lines with a "VERSION / TAG STATUS" header

### Requirement: Detect HEAD-already-tagged state
The hook SHALL detect when HEAD already has a version tag and display a "TAG OK" confirmation.

#### Scenario: HEAD is tagged with a version tag
- **WHEN** HEAD has one or more tags (e.g., `v0.3.0`)
- **THEN** the hook prints "TAG OK" followed by the tag name(s) and no further action is recommended

### Requirement: Detect VERSION-bumped-but-untagged state
The hook SHALL detect when the VERSION file contains a version that does not match the latest tag and HEAD is untagged, indicating the version was bumped but not yet tagged.

#### Scenario: VERSION file is ahead of latest tag
- **WHEN** HEAD is untagged AND the VERSION file value (e.g., `0.4.0`) differs from the latest tag's version (e.g., `v0.3.0` → `0.3.0`)
- **THEN** the hook prints "READY TO TAG" with the exact `git tag -a vX.Y.Z -m "vX.Y.Z"` command to run

### Requirement: Detect VERSION-not-bumped state
The hook SHALL detect when the VERSION file matches the latest tag but HEAD is untagged, indicating new commits exist without a version bump.

#### Scenario: VERSION matches latest tag and commits exist since tag
- **WHEN** HEAD is untagged AND the VERSION file value matches the latest tag's version AND there are commits since the latest tag
- **THEN** the hook prints "VERSION NOT BUMPED" with:
  - The number of commits since the last tag
  - Up to 10 commit subjects (one per line) since the last tag
  - If more than 10 commits exist, a note indicating how many total
  - Guidance to update the VERSION file and then tag

### Requirement: Show current state values
The hook SHALL always display the latest tag, the VERSION file value, and the commit count since last tag within the diagnostic section, so the reader has all facts at a glance.

#### Scenario: Diagnostic section includes state values
- **WHEN** the VERSION / TAG STATUS section is displayed
- **THEN** it includes labeled lines for: latest tag, VERSION file value, and number of commits since the latest tag

### Requirement: Graceful fallback when no tags exist
The hook SHALL skip the diagnostic section entirely when no tags exist in the repository.

#### Scenario: Repository has no tags
- **WHEN** `git describe --tags --abbrev=0` fails (no tags)
- **THEN** no VERSION / TAG STATUS section is printed
