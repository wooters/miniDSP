## ADDED Requirements

### Requirement: C17 standard compliance
All library source files, public headers, test files, and example programs SHALL compile cleanly under `-std=c17 -Wall -Wextra -pedantic` with zero warnings.

#### Scenario: Clean build with C17
- **WHEN** running `make clean && make` with `CFLAGS` containing `-std=c17`
- **THEN** the build completes with zero errors and zero warnings

#### Scenario: Test suite passes under C17
- **WHEN** running `make -C tests && ./tests/test_minidsp` after building with `-std=c17`
- **THEN** all existing tests pass with identical results

### Requirement: No C23 language features
The codebase SHALL NOT use any C23-only language features including `nullptr`, `[[attribute]]` syntax, or keyword `bool`/`true`/`false` without `<stdbool.h>`.

#### Scenario: No nullptr usage
- **WHEN** searching all `.c` and `.h` files for `nullptr`
- **THEN** zero matches are found

#### Scenario: stdbool.h included where bool is used
- **WHEN** any source file uses `bool`, `true`, or `false`
- **THEN** `<stdbool.h>` is included (directly or transitively via `minidsp.h`)

#### Scenario: No C23 attribute syntax
- **WHEN** searching all `.c` and `.h` files for `[[deprecated`
- **THEN** zero matches are found; `__attribute__((deprecated))` is used instead

### Requirement: Build configuration updated
The `config.mk` file SHALL specify `-std=c17` in `CFLAGS`.

#### Scenario: config.mk uses C17
- **WHEN** reading `config.mk`
- **THEN** `CFLAGS` contains `-std=c17` (not `-std=c23` or `-std=c2x`)
