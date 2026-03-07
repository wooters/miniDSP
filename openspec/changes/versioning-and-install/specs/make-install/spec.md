## ADDED Requirements

### Requirement: Install target copies library and headers
`make install` SHALL copy `libminidsp.a` to `$(PREFIX)/lib/` and all public headers from `include/` to `$(PREFIX)/include/`.

#### Scenario: Default prefix install
- **WHEN** running `make install` without setting PREFIX
- **THEN** files are installed to `/usr/local/lib/` and `/usr/local/include/`

#### Scenario: Custom prefix install
- **WHEN** running `make install PREFIX=/tmp/test-prefix`
- **THEN** `libminidsp.a` exists at `/tmp/test-prefix/lib/libminidsp.a`
- **THEN** all public headers exist under `/tmp/test-prefix/include/`

### Requirement: Uninstall target removes installed files
`make uninstall` SHALL remove the files installed by `make install` at the same `PREFIX`.

#### Scenario: Clean uninstall
- **WHEN** running `make uninstall PREFIX=/tmp/test-prefix` after a prior install to that prefix
- **THEN** `libminidsp.a` and all public headers are removed from that prefix

### Requirement: Install requires built library
`make install` SHALL depend on `libminidsp.a` so that the library is built before installation if needed.

#### Scenario: Install triggers build
- **WHEN** running `make install` on a clean checkout
- **THEN** the library is built first, then installed
