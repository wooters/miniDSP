## Context

miniDSP is a C17 library with 65 public API functions across 9 source files. It currently enforces preconditions exclusively via `assert()`, which calls `abort()` on failure — killing the host process. Assertions also vanish under `-DNDEBUG`, leaving no runtime protection in release builds.

The library has existing global mutable state (FFT plan caches, mel filterbank caches) managed through static variables, so a global error handler fits the existing architecture.

## Goals / Non-Goals

**Goals:**
- Every precondition violation reports an error and returns safely — never crashes the caller
- Checks are always active, even in release builds (`-DNDEBUG`)
- Callers can install a custom error handler for integration with their own logging/error systems
- Existing API signatures are preserved — no breaking changes to function return types or parameters
- All tools and examples demonstrate error handler installation

**Non-Goals:**
- Thread-safe error handler (documented as "set once before use", consistent with rest of library)
- Per-call error context or error stacks
- Changing sentinel return values for valid-but-unresolved runtime outcomes (e.g., `MD_f0_autocorrelation` returning 0.0 for "no peak found")
- Making the library broadly thread-safe

## Decisions

### 1. Global error handler callback (not return codes)

**Choice**: Configurable `MD_ErrorHandler` callback, not `int` return codes on every function.

**Alternatives considered**:
- **Return codes everywhere**: Would require changing all 39 `void` functions and 12 `double`-returning functions to return `int`, with output through pointers. Breaks the entire API, loses composability of `double`-returning functions, and callers can still ignore return values.
- **Thread-local handler**: `_Thread_local` gives per-thread handlers but has surprising semantics (main thread handler doesn't propagate) and portability issues (MSVC `__declspec(thread)`).

**Rationale**: Preserves all existing function signatures. This is the same pattern FFTW uses — miniDSP's primary dependency — so it's familiar to the target audience.

### 2. Internal check macros using `__func__`

**Choice**: `MD_CHECK(cond, code, msg, retval)` and `MD_CHECK_VOID(cond, code, msg)` macros in `minidsp_internal.h` that use `__func__` for automatic function name capture.

**Rationale**: `__func__` is guaranteed in C99+ (well within C17). Eliminates a manual string argument at every call site, reducing the chance of copy-paste errors where the function name string doesn't match the actual function.

### 3. Default handler logs to stderr

**Choice**: The default `MD_ErrorHandler` calls `fprintf(stderr, ...)`. Never aborts. Users who want strict/abort behavior can install a handler that calls `abort()`.

**Alternatives considered**:
- **Silent default**: Zero side effects but makes debugging hard — bad input silently produces wrong results.
- **Abort default (like FFTW)**: Defeats the purpose of replacing `assert()`.

**Rationale**: Stderr logging is the least surprising default. It's visible during development, and production callers can override or silence it.

### 4. New source file `src/minidsp_error.c`

**Choice**: Error handler state (`current_handler` pointer) and `md_report_error()` live in a new `src/minidsp_error.c`, not in an existing file.

**Rationale**: Error handling is orthogonal to all existing modules. Putting it in `minidsp_core.c` would conflate concerns. A dedicated file keeps the module boundaries clean.

### 5. Error codes as an enum, not `#define` constants

**Choice**: `typedef enum { MD_ERR_NULL_POINTER = 1, ... } MD_ErrorCode;`

**Rationale**: Enum gives type safety in the handler signature and debugger visibility. Starting at 1 reserves 0 for "no error" if ever needed.

### 6. Safe default return values

**Choice**: On precondition failure, functions return type-appropriate safe defaults:
- `void` → early return (no-op)
- `double` → `0.0`
- `unsigned` → `0`
- `int` → `-1` (matches existing `fileio.c` convention)

**Rationale**: These defaults are the least harmful values. A caller who doesn't check for errors gets a zero/no-op rather than a crash or garbage data.

## Risks / Trade-offs

**[Silent wrong results]** → A void function that no-ops on bad input leaves the caller's output buffer unchanged. The caller may not realize the computation didn't happen. **Mitigation**: Default handler logs to stderr, making failures visible. Documentation emphasizes checking handler output during development.

**[Performance overhead]** → Every function call now has `if` checks instead of assertions that compile away. **Mitigation**: Branch prediction makes not-taken `if` checks essentially free (~1 cycle). For DSP inner loops processing thousands of samples, the overhead of a handful of pointer/size checks at function entry is negligible.

**[Macro hygiene]** → `MD_CHECK` macros use `do { ... } while(0)` and evaluate arguments once. `__func__` is a compiler-provided identifier, not a macro, so there are no expansion issues. **Risk is low**.

**[Global mutable state]** → The handler pointer is global, same as FFT plan caches. **Mitigation**: Documented as "set once before use". Upgrading to `_Atomic` later is a one-line change if thread safety becomes a goal.
