## 1. Version Infrastructure

- [x] 1.1 Create `VERSION` file at repo root containing `0.1.0`
- [x] 1.2 Add version macros to `include/minidsp.h` with `#ifndef` guards for default values
- [x] 1.3 Update root `Makefile` to read `VERSION` and pass `-D` flags via CFLAGS

## 2. Install/Uninstall Targets

- [x] 2.1 Add `install` target to root Makefile (copies `libminidsp.a` and `include/*.h` to PREFIX)
- [x] 2.2 Add `uninstall` target to root Makefile
- [x] 2.3 Add `PREFIX` variable with `/usr/local` default

## 3. Verification

- [x] 3.1 Run `make clean && make` — verify version macros compile without warnings
- [x] 3.2 Run `make install PREFIX=/tmp/minidsp-test` — verify files are installed correctly
- [x] 3.3 Run `make uninstall PREFIX=/tmp/minidsp-test` — verify files are removed
- [x] 3.4 Run tests to confirm nothing is broken

## 4. Documentation and Tagging

- [x] 4.1 Update `CLAUDE.md` with versioning conventions
- [ ] 4.2 Create annotated git tag `v0.1.0` (after commit)
