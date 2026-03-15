## ADDED Requirements

### Requirement: Shared test infrastructure header
The test suite SHALL provide a `tests/test_helpers.h` header that contains `extern` declarations for the global test counters (`tests_run`, `tests_passed`, `tests_failed`), the `RUN_TEST` macro, the `approx_equal` static inline function, and the `delay_signal` static helper.

#### Scenario: Header is self-contained
- **WHEN** any `test_*.c` file includes only `test_helpers.h` and the library headers it needs
- **THEN** it SHALL compile without errors or missing symbol warnings

#### Scenario: Counters are shared across files
- **WHEN** multiple test files each call `RUN_TEST` and the driver prints the summary
- **THEN** the reported totals SHALL reflect the sum of all tests across all files

### Requirement: Per-module test files
Each source module SHALL have a corresponding test file in `tests/` containing all tests for that module's API functions. Test functions within each file SHALL remain `static`. Each file SHALL expose exactly one non-static function with the signature `void run_<module>_tests(void)`.

#### Scenario: Test file exists for every module
- **WHEN** a developer looks for tests for `MD_foo` defined in `src/minidsp_xxx.c`
- **THEN** the tests SHALL be in `tests/test_xxx.c`

#### Scenario: No test function name collisions
- **WHEN** all test files are compiled and linked together
- **THEN** there SHALL be no duplicate symbol errors (test functions are `static`)

### Requirement: Driver file orchestration
`tests/test_minidsp.c` SHALL be a slim driver that defines the global test counters, calls each module's `run_*_tests()` function, calls `MD_shutdown()`, and prints the pass/fail summary. It SHALL NOT contain any test logic.

#### Scenario: Driver runs all tests
- **WHEN** `./test_minidsp` is executed
- **THEN** it SHALL run all 277 tests and print the same summary format: `=== Results: N/N passed ===`

#### Scenario: Driver exit code
- **WHEN** any test fails
- **THEN** the driver SHALL exit with code 1
- **WHEN** all tests pass
- **THEN** the driver SHALL exit with code 0

### Requirement: Build system compiles all test files
`tests/Makefile` SHALL compile each `test_*.c` file to a `.o` and link them all into the `test_minidsp` binary. The `make test` target SHALL continue to build and run the suite.

#### Scenario: Incremental rebuild
- **WHEN** a developer modifies only `tests/test_generators.c`
- **THEN** `make test_minidsp` SHALL recompile only `test_generators.o` and re-link

#### Scenario: Clean build
- **WHEN** `make clean` is run in `tests/`
- **THEN** all `.o` files and the `test_minidsp` binary SHALL be removed

### Requirement: Zero behavior change
The refactored test suite SHALL produce identical test results — same tests, same names, same output format, same pass/fail outcomes. No tests SHALL be added, removed, or modified in logic.

#### Scenario: Full pass count preserved
- **WHEN** `make test` runs before and after the refactor
- **THEN** the number of tests run and passed SHALL be identical (277/277)

#### Scenario: Output format preserved
- **WHEN** the test suite runs
- **THEN** section headers (`--- MD_xxx ---`) and per-test `[PASS]`/`[FAIL]` lines SHALL appear in the same order as before
