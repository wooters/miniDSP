## Context

The miniDSP library currently compiles with `-std=c23`. Three C23 features are used:

1. **`nullptr`** — 191 occurrences across 15 files (src/, include/, tests/, examples/)
2. **`bool`/`true`/`false` as keywords** — 3 occurrences in `minidsp.h` and `minidsp_core.c` (no `<stdbool.h>` include)
3. **`[[deprecated("msg")]]`** — 1 occurrence in `fileio.h`

No other C23 features (`constexpr`, `typeof`, `_BitInt`, `#embed`, `auto` type inference, `static_assert` without message) are used.

## Goals / Non-Goals

**Goals:**
- Compile cleanly with `-std=c17 -Wall -Wextra -pedantic` on GCC 9+ and Clang 10+
- Mechanical, behavior-preserving replacements only
- Update project documentation to reflect the new C17 baseline

**Non-Goals:**
- Removing `-pedantic` (keep strict compliance checking)
- Supporting C11 or older standards
- Changing any runtime behavior or public API signatures

## Decisions

### 1. `nullptr` → `NULL`

**Decision**: Replace all `nullptr` with `NULL`.

**Rationale**: `NULL` is the universal C null pointer constant, available in `<stddef.h>` (already included transitively via other standard headers in `minidsp.h`). No semantic difference for pointer comparisons and assignments.

**Alternative considered**: Define `#define nullptr NULL` in a compat header — adds complexity for no benefit since the replacement is mechanical.

### 2. `bool` / `true` / `false` → `<stdbool.h>`

**Decision**: Add `#include <stdbool.h>` to `minidsp.h` (which all source files include).

**Rationale**: In C17, `bool`, `true`, `false` are macros from `<stdbool.h>`, not keywords. Adding the include makes them available everywhere with zero code changes to usage sites.

### 3. `[[deprecated]]` → `__attribute__((deprecated))`

**Decision**: Replace `[[deprecated("msg")]]` with `__attribute__((deprecated("msg")))`.

**Rationale**: The GCC/Clang attribute syntax is supported by all target compilers (GCC 4+, Clang 3+) and is the standard pre-C23 approach. Only one occurrence exists.

### 4. Build flag change

**Decision**: Change `CFLAGS` in `config.mk` from `-std=c23` to `-std=c17`.

**Rationale**: Single point of change — all three Makefiles include `config.mk`.

## Risks / Trade-offs

- **[Low] Future C23 feature adoption blocked** → Acceptable trade-off for CI compatibility. Can revisit when C23 support is widespread in CI environments.
- **[Low] Missed `nullptr` occurrence** → Mitigated by grep verification after replacement and a clean compile test.
- **[Low] `NULL` in pointer arithmetic contexts** → Not a concern; all current `nullptr` usage is assignment/comparison, never arithmetic.
