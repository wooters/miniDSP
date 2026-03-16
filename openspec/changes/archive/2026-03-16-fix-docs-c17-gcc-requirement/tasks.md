## 1. README Updates

- [x] 1.1 Change the language badge from "C23" to "C17" and update the Wikipedia link (line 6)
- [x] 1.2 Change `-std=c23` to `-std=c17` in the direct compile command (line 146)
- [x] 1.3 Change `-std=c23` to `-std=c17` in the Makefile example (line 156)
- [x] 1.4 Remove the paragraph about needing GCC 14 on Ubuntu (line 439)
- [x] 1.5 Update the container-test description to remove "with GCC 14" (line 455)

## 2. Dockerfile Simplification

- [x] 2.1 Replace `gcc-14` with `gcc` in the apt-get install line
- [x] 2.2 Remove the `update-alternatives` command

## 3. Verification

- [x] 3.1 Run `make` to confirm the library still builds
- [x] 3.2 Run `make container-test` to confirm the container builds and tests pass with default GCC
