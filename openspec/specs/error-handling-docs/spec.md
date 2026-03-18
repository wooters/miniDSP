## ADDED Requirements

### Requirement: Error handling API documentation
The public API additions (`MD_ErrorCode`, `MD_ErrorHandler`, `MD_set_error_handler()`) SHALL have Doxygen doc-comments in `minidsp.h` that include:
- Description of each error code
- The handler function signature and when it is called
- Thread-safety contract ("set once before use")
- Example code showing handler installation

#### Scenario: Doxygen generates error handling reference
- **WHEN** `make docs` is run
- **THEN** the generated HTML SHALL include documentation for `MD_ErrorCode`, `MD_ErrorHandler`, and `MD_set_error_handler()` with descriptions and usage examples

### Requirement: Error contract documented in library overview
The library's Doxygen main page or a dedicated guide page SHALL document the error handling contract:
- The library never aborts on precondition violations
- Default behavior: log to stderr and return safe defaults
- How to install a custom handler
- List of error codes and when each is used
- Thread-safety constraints

#### Scenario: User finds error handling guide
- **WHEN** a user browses the miniDSP Doxygen documentation
- **THEN** they SHALL find a section or page explaining the error handling contract with code examples

### Requirement: CLAUDE.md updated
The API design contract in `CLAUDE.md` SHALL be updated to reflect:
- `assert()` replaced by `MD_CHECK` / `MD_CHECK_VOID` macros
- Error handler mechanism and default behavior
- Safe default return value convention

#### Scenario: CLAUDE.md reflects new contract
- **WHEN** a developer reads the CLAUDE.md API design contract section
- **THEN** it SHALL describe the error handler system, not `assert()`
