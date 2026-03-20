## 1. Doxygen Configuration

- [x] 1.1 Add `GENERATE_XML = YES` to `Doxyfile`
- [x] 1.2 Update `make clean` in root Makefile to also remove `docs/xml/`

## 2. Generation Script

- [x] 2.1 Create `scripts/gen_llms_txt.py` — parse Doxygen XML (`docs/xml/minidsp_8h.xml`) to extract public function signatures, descriptions, parameters, return values, and code examples
- [x] 2.2 Add guide markdown processing — read `guides/*.md` files, strip Doxygen syntax (`\f$`→`$`, remove `\htmlonly` blocks, strip `{#id}` anchors, convert `\subpage`/`\ref` to plain text)
- [x] 2.3 Add `\snippet` resolution — locate source files in `EXAMPLE_PATH` directories, extract code between `//! [id]` markers, inline as fenced code blocks
- [x] 2.4 Generate `llms.txt` index file (library name, description, links to llms-full.txt and repo)
- [x] 2.5 Generate `llms-full.txt` combining API reference + guide content, with guides ordered by `\subpage` directives in `guides/tutorials.md`
- [x] 2.6 Add error handling: exit non-zero with message if `docs/xml/` directory is missing

## 3. Build Pipeline Integration

- [x] 3.1 Add `python3 scripts/gen_llms_txt.py` to the `docs` target in root Makefile (after the `doxygen Doxyfile` line)
- [x] 3.2 Add `<link rel="help" type="text/markdown" href="llms-full.txt">` to `doxygen-custom/header.html` `<head>` section

## 4. Verification

- [x] 4.1 Run `make docs` end-to-end and verify `docs/html/llms.txt` and `docs/html/llms-full.txt` are generated
- [x] 4.2 Verify `llms-full.txt` contains all public API functions from `minidsp.h`
- [x] 4.3 Verify guide content is included with Doxygen syntax properly stripped
- [x] 4.4 Verify `make clean` removes `docs/xml/`
