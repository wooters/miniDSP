## Why

miniDSP currently uses `assert()` for precondition checking (244 assertions across 9 source files). This is unsafe for a library: `assert()` calls `abort()`, killing the caller's entire process on bad input. Worse, assertions vanish under `-DNDEBUG`, leaving undefined behavior in release builds. A library should never crash its host application.

## What Changes

- **New error handler system**: Configurable callback that fires on precondition violations, with a default handler that logs to stderr and returns safely — never aborts
- **BREAKING**: All 244 `assert()` calls replaced with always-on `MD_CHECK` / `MD_CHECK_VOID` macros that report errors and return safe defaults (0.0, 0, no-op) instead of aborting
- **New public API**: `MD_ErrorCode` enum, `MD_ErrorHandler` typedef, `MD_set_error_handler()` function
- **Tools and examples updated**: All programs in `tools/` and `examples/` install an appropriate error handler
- **Documentation**: Doxygen docs for the error handling API and the library's error contract
- **Tests**: Test coverage for handler installation, error reporting, and safe default returns

## Capabilities

### New Capabilities
- `error-handling`: Configurable error handler callback system — error codes, check macros, default stderr handler, safe early returns on precondition failure
- `error-handling-docs`: Doxygen documentation for the error handling API and the library's error contract

### Modified Capabilities

## Impact

- **Public API**: New types (`MD_ErrorCode`, `MD_ErrorHandler`) and function (`MD_set_error_handler`) added to `minidsp.h`
- **All source files**: 9 `.c` files modified to replace `assert()` with `MD_CHECK` macros
- **Internal header**: `minidsp_internal.h` gains `MD_CHECK` / `MD_CHECK_VOID` macro definitions and `md_report_error()` declaration
- **New source file**: `src/minidsp_error.c` for handler state and `md_report_error()` implementation
- **Build system**: Root Makefile updated to compile `minidsp_error.c` into `libminidsp.a`
- **Tools and examples**: All programs updated to call `MD_set_error_handler()`
- **Tests**: New test functions in `test_minidsp.c` for error handling behavior
- **CLAUDE.md**: API design contract updated to reflect new error handling approach
