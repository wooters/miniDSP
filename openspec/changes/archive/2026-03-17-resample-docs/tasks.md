## 1. Guide Page

- [x] 1.1 Create `guides/resample-tool.md` with anchor `{#resample-tool}`: overview, build/run instructions, options table (`-z`, `-b` with defaults and descriptions), example workflows, link to `\ref resampling` for math detail
- [x] 1.2 Add snippet markers (`//! [id]`) in `tools/resample/resample.c` for any code to embed via `\snippet`

## 2. Doxygen Integration

- [x] 2.1 Add `\subpage resample-tool` to `guides/tutorials.md` under the Tools section
- [x] 2.2 Add `tools/resample` to `EXAMPLE_PATH` in `Doxyfile`

## 3. Verification

- [x] 3.1 Run `make docs` and verify the resample tool guide renders correctly in `docs/html/`
