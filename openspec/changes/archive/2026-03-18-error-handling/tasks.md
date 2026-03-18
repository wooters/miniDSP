## 1. Error Handler Infrastructure

- [x] 1.1 Add `MD_ErrorCode` enum and `MD_ErrorHandler` typedef to `include/minidsp.h` with Doxygen doc-comments
- [x] 1.2 Add `MD_set_error_handler()` declaration to `include/minidsp.h` with Doxygen doc-comment (including thread-safety contract and usage example)
- [x] 1.3 Create `src/minidsp_error.c` with default handler, `current_handler` static, `MD_set_error_handler()` implementation, and `md_report_error()` implementation
- [x] 1.4 Add `MD_CHECK` and `MD_CHECK_VOID` macro definitions to `src/minidsp_internal.h`, plus `md_report_error()` declaration
- [x] 1.5 Update root `Makefile` to compile `minidsp_error.c` into `libminidsp.a`

## 2. Replace assert() in Source Files

- [x] 2.1 Replace all `assert()` calls in `src/minidsp_core.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.2 Replace all `assert()` calls in `src/minidsp_spectrum.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.3 Replace all `assert()` calls in `src/minidsp_fir.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.4 Replace all `assert()` calls in `src/minidsp_generators.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.5 Replace all `assert()` calls in `src/minidsp_steg.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.6 Replace all `assert()` calls in `src/minidsp_dtmf.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.7 Replace all `assert()` calls in `src/minidsp_resample.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.8 Replace all `assert()` calls in `src/minidsp_spectext.c` with `MD_CHECK` / `MD_CHECK_VOID` macros
- [x] 2.9 Replace all `assert()` calls in `src/minidsp_gcc.c` with `MD_CHECK` / `MD_CHECK_VOID` macros

## 3. Tests

- [x] 3.1 Add error handling tests to `tests/test_minidsp.c`: custom handler installation, default handler restoration, correct error code/function name/message delivery
- [x] 3.2 Add safe-default-return tests: verify `double` functions return `0.0`, `void` functions no-op, `unsigned` functions return `0` on bad input
- [x] 3.3 Verify all existing tests still pass after the `assert()` → `MD_CHECK` migration

## 4. Tools and Examples

- [x] 4.1 Update all example programs in `examples/` to call `MD_set_error_handler()` at startup
- [x] 4.2 Update all tool programs in `tools/` to call `MD_set_error_handler()` at startup

## 5. Documentation

- [x] 5.1 Add or update Doxygen guide page documenting the error handling contract (error codes, default behavior, custom handler installation, thread-safety, code examples)
- [x] 5.2 Update the `@brief` feature list in `include/minidsp.h` to mention configurable error handling
- [x] 5.3 Update `CLAUDE.md` API design contract to reflect `MD_CHECK` macros and error handler system (replacing the `assert()` convention)
