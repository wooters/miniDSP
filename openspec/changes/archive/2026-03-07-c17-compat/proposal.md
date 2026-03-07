## Why

The codebase requires `-std=c23`, which is only supported by very recent compilers (GCC 14+, Clang 17+). Many CI pipelines (e.g., Ubuntu 22.04 runners, older Docker images) ship older toolchains that lack C23 support. Downgrading to C17 broadens compatibility without losing meaningful functionality.

## What Changes

- **BREAKING**: Change `-std=c23` to `-std=c17` in `config.mk`
- Replace all `nullptr` usage (191 occurrences, 15 files) with `NULL`
- Add `#include <stdbool.h>` where `bool`/`true`/`false` are used (C17 requires this header)
- Replace `[[deprecated("msg")]]` C23 attribute with `__attribute__((deprecated("msg")))` (GCC/Clang compatible)
- Update `CLAUDE.md` C23 notes to reflect the new C17 baseline

## Capabilities

### New Capabilities
- `c17-compat`: Ensure all library source, headers, tests, and examples compile cleanly under `-std=c17 -Wall -Wextra -pedantic`

### Modified Capabilities
<!-- No existing specs to modify -->

## Impact

- **Build system**: `config.mk` CFLAGS change from `-std=c23` to `-std=c17`
- **All source files**: Mechanical `nullptr` → `NULL` replacement across `src/`, `include/`, `tests/`, `examples/`
- **Headers**: `<stdbool.h>` include added to `minidsp.h`
- **Public API**: `[[deprecated]]` attribute syntax changes (no functional change)
- **Documentation**: `CLAUDE.md` C23-specific notes updated
- **CI compatibility**: Enables builds on GCC 9+, Clang 10+, and standard Ubuntu 20.04+ runners
