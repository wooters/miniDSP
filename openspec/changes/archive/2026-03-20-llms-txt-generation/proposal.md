## Why

AI coding agents (Claude Code, Cursor, Copilot, etc.) increasingly need to understand library APIs to generate correct code. Doxygen HTML is terrible for this — fragmented across dozens of pages with navigation chrome, JavaScript, and CSS. The emerging `llms.txt` convention (analogous to `robots.txt`) solves this by serving clean, structured markdown at well-known paths that agents can ingest in a single `fetch` or `curl` call.

miniDSP's docs are already comprehensive (21 guides, full API doc-comments with formulas and examples). We just need a build step to extract that content into agent-friendly plain text and wire it into the existing deployment pipeline so it stays current automatically.

## What Changes

- Enable Doxygen XML output (`GENERATE_XML = YES`) as structured source data for the generation script
- Add a Python script (`scripts/gen_llms_txt.py`) that parses Doxygen XML and guide markdown to produce `llms.txt` (concise index) and `llms-full.txt` (complete API reference + guides)
- Extend the `make docs` target to run the generation script after Doxygen, placing output in `docs/html/` for automatic GitHub Pages deployment
- Install Python 3 in the CI workflow (`.github/workflows/docs.yml`) and run the script as part of the docs build
- Add a `<link>` tag in `doxygen-custom/header.html` for programmatic discoverability

## Capabilities

### New Capabilities
- `llms-txt-gen`: Generation of `llms.txt` and `llms-full.txt` from Doxygen XML output and guide markdown files, integrated into the docs build pipeline

### Modified Capabilities
<!-- None — this is purely additive -->

## Impact

- **Doxyfile**: One-line addition (`GENERATE_XML = YES`)
- **Root Makefile**: `docs` target gains a Python script invocation after `doxygen`
- **CI workflow**: Needs `python3` available (already present on `ubuntu-24.04`) and the script run
- **doxygen-custom/header.html**: Optional `<link>` tag for discoverability
- **New files**: `scripts/gen_llms_txt.py`, output `docs/html/llms.txt` and `docs/html/llms-full.txt` (generated, not tracked)
- **Dependencies**: Python 3 standard library only (xml.etree.ElementTree for XML parsing) — no pip packages
