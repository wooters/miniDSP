## ADDED Requirements

### Requirement: Error code enumeration
The library SHALL define an `MD_ErrorCode` enum with the following codes:
- `MD_ERR_NULL_POINTER = 1` — a required pointer argument is NULL
- `MD_ERR_INVALID_SIZE` — a size or count argument is invalid (e.g., N == 0)
- `MD_ERR_INVALID_RANGE` — a range or bound argument is invalid (e.g., min >= max, frequency out of bounds)
- `MD_ERR_ALLOC_FAILED` — a memory allocation failed

The enum SHALL start at 1, reserving 0 for "no error".

#### Scenario: Error codes are distinct nonzero values
- **WHEN** any error code from `MD_ErrorCode` is examined
- **THEN** it SHALL be a positive integer, distinct from all other error codes

### Requirement: Error handler callback type
The library SHALL define a function pointer type `MD_ErrorHandler` with signature:
```c
typedef void (*MD_ErrorHandler)(MD_ErrorCode code, const char *func_name, const char *message);
```

#### Scenario: Handler receives error details
- **WHEN** a precondition violation occurs in `MD_energy()` with a NULL signal pointer
- **THEN** the handler SHALL be called with `code == MD_ERR_NULL_POINTER`, `func_name` equal to `"MD_energy"`, and `message` describing the violation

### Requirement: Custom error handler installation
The library SHALL provide `void MD_set_error_handler(MD_ErrorHandler handler)` that installs a custom error handler. Passing NULL SHALL restore the default handler.

#### Scenario: Install custom handler
- **WHEN** a caller installs a custom handler via `MD_set_error_handler(my_handler)`
- **THEN** all subsequent precondition violations SHALL invoke `my_handler`

#### Scenario: Restore default handler
- **WHEN** a caller calls `MD_set_error_handler(NULL)`
- **THEN** the default stderr-logging handler SHALL be restored

#### Scenario: Set-once-before-use threading contract
- **WHEN** `MD_set_error_handler()` is called
- **THEN** it MUST be called before any other `MD_*` function and MUST NOT be called concurrently with any other `MD_*` function

### Requirement: Default error handler
The library SHALL provide a default error handler that logs the error to stderr in the format:
```
miniDSP: error <code> in <func_name>(): <message>
```
The default handler SHALL NOT call `abort()` or terminate the process.

#### Scenario: Default handler logs to stderr
- **WHEN** a precondition violation occurs and no custom handler has been installed
- **THEN** a human-readable error message SHALL be written to stderr

#### Scenario: Default handler does not abort
- **WHEN** a precondition violation occurs with the default handler active
- **THEN** the calling process SHALL NOT be terminated

### Requirement: Check macros replace assert
All `assert()` calls in the library SHALL be replaced with `MD_CHECK` or `MD_CHECK_VOID` macros that:
1. Evaluate the condition
2. On failure: call the error handler via `md_report_error()`, then return a safe default
3. Are always active — they MUST NOT be disabled by `-DNDEBUG`

#### Scenario: Check fires in release build
- **WHEN** a library is compiled with `-DNDEBUG` and a precondition violation occurs
- **THEN** the error handler SHALL still be called and the function SHALL return a safe default

#### Scenario: Check passes on valid input
- **WHEN** all preconditions are satisfied
- **THEN** the check macros SHALL have no observable side effects and the function SHALL execute normally

### Requirement: Safe default return values
On precondition failure, functions SHALL return type-appropriate safe defaults:
- `void` functions: early return (no-op)
- `double`-returning functions: return `0.0`
- `unsigned`-returning functions: return `0`
- `int`-returning functions: return `-1`

#### Scenario: Void function returns safely on bad input
- **WHEN** `MD_hanning_window(NULL, 512)` is called
- **THEN** the error handler SHALL fire and the function SHALL return without writing to any buffer

#### Scenario: Double function returns 0.0 on bad input
- **WHEN** `MD_energy(NULL, 512)` is called
- **THEN** the error handler SHALL fire and the function SHALL return `0.0`

#### Scenario: Unsigned function returns 0 on bad input
- **WHEN** `MD_stft_num_frames(NULL, 0, 0, 0)` is called with invalid parameters
- **THEN** the error handler SHALL fire and the function SHALL return `0`

### Requirement: Sentinel returns unchanged
Existing sentinel return values for valid-but-unresolved runtime outcomes (e.g., `MD_f0_autocorrelation` returning `0.0` for "no peak found") SHALL remain unchanged. These are not errors and MUST NOT trigger the error handler.

#### Scenario: No-peak-found is not an error
- **WHEN** `MD_f0_autocorrelation()` is called with valid input but no autocorrelation peak is found
- **THEN** it SHALL return `0.0` without invoking the error handler

### Requirement: All public functions checked
Every public `MD_*` function SHALL validate its preconditions using the check macros. At minimum, every function SHALL check:
- All pointer arguments for NULL
- All size/count arguments for validity (e.g., N > 0)
- All range arguments for consistency (e.g., min < max)

#### Scenario: NULL pointer caught in every function
- **WHEN** any public `MD_*` function that takes a pointer argument is called with NULL for that argument
- **THEN** the error handler SHALL be invoked with `MD_ERR_NULL_POINTER`

### Requirement: Tools and examples install error handler
All programs in `tools/` and `examples/` SHALL call `MD_set_error_handler()` at startup to demonstrate proper usage.

#### Scenario: Example program installs handler
- **WHEN** any example program in `examples/` starts execution
- **THEN** it SHALL have called `MD_set_error_handler()` before any other `MD_*` calls

#### Scenario: Tool program installs handler
- **WHEN** any tool program in `tools/` starts execution
- **THEN** it SHALL have called `MD_set_error_handler()` before any other `MD_*` calls

### Requirement: Error handling tests
The test suite SHALL include tests that verify:
- Custom handler installation and invocation
- Default handler restoration via `MD_set_error_handler(NULL)`
- Safe default return values for each return type category
- Error handler receiving correct error code, function name, and message

#### Scenario: Test verifies handler receives correct error code
- **WHEN** a test installs a custom handler and calls `MD_energy(NULL, 100)`
- **THEN** the custom handler SHALL have been called with `code == MD_ERR_NULL_POINTER`

#### Scenario: Test verifies safe return value
- **WHEN** a test calls `MD_energy(NULL, 100)`
- **THEN** the return value SHALL be `0.0`
