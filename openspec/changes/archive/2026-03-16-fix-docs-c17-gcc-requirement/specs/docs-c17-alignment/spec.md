## ADDED Requirements

### Requirement: README reflects C17 standard
All references to the C language standard in README.md SHALL state C17, matching the actual `-std=c17` flag in `config.mk`.

#### Scenario: Language badge shows C17
- **WHEN** a user views the README badges
- **THEN** the language badge SHALL display "C17" and link to the C17 Wikipedia page

#### Scenario: Compile examples use C17 flag
- **WHEN** a user follows the "Use in your project" compile instructions
- **THEN** the `-std=` flag SHALL be `-std=c17` in both the direct command and the Makefile example

### Requirement: No GCC 14 requirement for Ubuntu
The README SHALL NOT state that GCC 14 is required on Ubuntu, since C17 is supported by the default system compiler (GCC 13) on Ubuntu 24.04.

#### Scenario: Ubuntu build instructions
- **WHEN** a user reads the "Build and Test" section on Ubuntu 24.04
- **THEN** there SHALL be no instruction to install `gcc-14` or any specific GCC version

### Requirement: Dockerfile uses system default GCC
The Dockerfile SHALL use the default `gcc` package instead of explicitly installing `gcc-14`.

#### Scenario: Container build with default GCC
- **WHEN** `make container-test` builds the Docker image
- **THEN** the image SHALL install `gcc` (not `gcc-14`) and SHALL NOT run `update-alternatives`
