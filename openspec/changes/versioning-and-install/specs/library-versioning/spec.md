## ADDED Requirements

### Requirement: VERSION file exists at repo root
The repository SHALL contain a `VERSION` file at the root containing a single line with the semver version string (e.g., `0.1.0`), with no leading `v` prefix.

#### Scenario: VERSION file format
- **WHEN** reading the `VERSION` file
- **THEN** it contains exactly one line matching the pattern `MAJOR.MINOR.PATCH` where each component is a non-negative integer

### Requirement: Compile-time version macros
The public header `minidsp.h` SHALL provide version macros: `MINIDSP_VERSION` (string), `MINIDSP_VERSION_MAJOR`, `MINIDSP_VERSION_MINOR`, `MINIDSP_VERSION_PATCH` (integers).

#### Scenario: Version macros match VERSION file
- **WHEN** building with the project Makefile and printing `MINIDSP_VERSION` from a C program
- **THEN** the printed string matches the content of the `VERSION` file

#### Scenario: Header compiles without Makefile
- **WHEN** compiling a C file that includes `minidsp.h` without the Makefile's `-D` flags
- **THEN** the version macros still have default values and compilation succeeds

### Requirement: Git tags follow semver
Release commits SHALL be tagged with annotated git tags in the format `vMAJOR.MINOR.PATCH`.

#### Scenario: Tag matches VERSION file
- **WHEN** checking out a release tag
- **THEN** the `VERSION` file content matches the tag name without the `v` prefix
