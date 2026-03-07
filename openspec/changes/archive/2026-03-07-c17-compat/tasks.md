## 1. Build Configuration

- [x] 1.1 Change `-std=c23` to `-std=c17` in `config.mk`

## 2. Replace C23 Features

- [x] 2.1 Add `#include <stdbool.h>` to `include/minidsp.h`
- [x] 2.2 Replace `[[deprecated("msg")]]` with `__attribute__((deprecated("msg")))` in `include/fileio.h`
- [x] 2.3 Replace all `nullptr` with `NULL` in `include/` headers
- [x] 2.4 Replace all `nullptr` with `NULL` in `src/` source files
- [x] 2.5 Replace all `nullptr` with `NULL` in `tests/test_minidsp.c`
- [x] 2.6 Replace all `nullptr` with `NULL` in `examples/` source files

## 3. Verification

- [x] 3.1 Grep for remaining `nullptr` — expect zero matches
- [x] 3.2 Grep for remaining `[[deprecated` — expect zero matches
- [x] 3.3 Run `make clean && make` — expect zero errors and zero warnings
- [x] 3.4 Run `make -C tests clean && make -C tests && ./tests/test_minidsp` — all tests pass
- [x] 3.5 Run `make -C examples clean && make -C examples` — all examples build

## 4. Documentation

- [x] 4.1 Update `CLAUDE.md` C23 notes to reflect C17 baseline
