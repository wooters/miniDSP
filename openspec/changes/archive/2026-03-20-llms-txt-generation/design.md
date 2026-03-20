## Context

miniDSP generates HTML documentation via Doxygen, deployed to GitHub Pages on every push to `main`. The Doxygen pipeline already generates audio samples, interactive plots, and HTML — but produces no plain-text output suitable for AI coding agents.

The public API is documented in `include/minidsp.h` with Doxygen doc-comments (formulas, `@param`, `@return`, `@code` examples). There are 21 guide markdown files in `guides/` using Doxygen-specific syntax (`\f$...\f$` for math, `\snippet` for code embedding, `\htmlonly` blocks for media).

The CI workflow (`.github/workflows/docs.yml`) runs `make docs` on Ubuntu 24.04 with Python 3 already available in the base image.

## Goals / Non-Goals

**Goals:**
- Generate `llms.txt` (concise library index) and `llms-full.txt` (complete API reference + guide content) as part of every docs build
- Output lands in `docs/html/` so GitHub Pages deploys it automatically — no separate workflow needed
- Use only Python 3 standard library (no pip dependencies)
- Strip Doxygen-specific syntax from guide content so output is clean markdown readable by any LLM

**Non-Goals:**
- Per-page markdown downloads or a markdown API server
- JSON or structured-data API reference format
- Versioned docs (multiple versions of llms.txt) — single latest version is sufficient
- Modifying any existing Doxygen HTML output or theme behavior

## Decisions

### 1. Doxygen XML as the API data source

**Decision**: Enable `GENERATE_XML = YES` and parse the XML output with Python's `xml.etree.ElementTree`.

**Why**: The XML output contains structured function signatures, parameter descriptions, return values, and detailed descriptions with formula markup — all machine-parseable. Parsing `minidsp.h` directly would require a C comment parser and miss Doxygen's cross-referencing.

**Alternative considered**: Regex parsing of `minidsp.h` doc-comments. Simpler but fragile — would break on multi-line comments, nested formulas, or `@code` blocks with special characters.

### 2. Guide markdown: direct file reading with syntax stripping

**Decision**: Read guide `.md` files directly from `guides/` and strip Doxygen syntax via regex transformations.

**Why**: The guide files are already markdown — the Doxygen XML representation of markdown pages loses formatting context. Direct file reading preserves the authored structure while a small set of regex transforms handles the Doxygen-specific syntax:
- `\f$...\f$` → `$...$` (inline math)
- `\f[...\f]` → `$$...$$` (display math)
- `\htmlonly`...`\endhtmlonly` blocks → removed (audio/iframe embeds aren't useful in text)
- `\snippet filename.c id` → resolved to inline code from the actual source file using `//! [id]` markers
- `{#anchor-id}` suffixes on headings → stripped
- `\subpage`, `\ref` directives → converted to plain text labels

### 3. Two-file output: llms.txt (index) + llms-full.txt (complete)

**Decision**: Generate two files following the llms.txt convention.

**Why**: `llms.txt` is small enough for agent discovery and routing (< 1KB). `llms-full.txt` contains everything an agent needs in one fetch — function signatures, descriptions, parameter docs, formulas, code examples, and tutorial content. An agent building with miniDSP fetches `llms-full.txt` once and has full context.

### 4. Script runs after Doxygen in the same `make docs` target

**Decision**: Add `python3 scripts/gen_llms_txt.py` to the `docs` target in the root Makefile, after the `doxygen Doxyfile` line.

**Why**: The script depends on Doxygen XML output (`docs/xml/`). Running it in the same target ensures correct ordering and means no CI workflow changes beyond what `make docs` already provides. Python 3 is pre-installed on `ubuntu-24.04`.

### 5. Snippet resolution from source files

**Decision**: The script resolves `\snippet filename.c snippet-id` directives by reading the referenced source file and extracting content between `//! [snippet-id]` markers.

**Why**: Snippets are a core part of the guide content — they show real, working code examples. Leaving them as `\snippet` references would be useless to an agent. The `EXAMPLE_PATH` directories in the Doxyfile (`examples`, `tools/audio_steg`, `tools/resample`) tell us where to find the source files.

### 6. Header link tag for discoverability

**Decision**: Add `<link rel="help" type="text/markdown" href="llms-full.txt">` to `doxygen-custom/header.html`.

**Why**: Agents that crawl a docs site can discover the llms.txt files programmatically via the `<link>` tag, similar to how RSS feeds are discovered. Low-effort addition with no visual impact.

## Risks / Trade-offs

- **Doxygen XML format changes**: Doxygen XML schema isn't formally versioned. Risk is low since we pin Doxygen 1.16.1 in CI → Mitigation: the script targets basic elements (`memberdef`, `param`, `detaileddescription`) that have been stable across Doxygen versions.

- **Large output file**: `llms-full.txt` could grow as guides are added. Currently estimated at ~50-80KB (21 guides + full API reference) → Mitigation: well within LLM context windows (100K+ tokens). If it grows past 200KB in the future, we can add section splitting.

- **Snippet resolution fragility**: If snippet markers are renamed or files moved, the script will silently omit the code block → Mitigation: the script prints a warning to stderr for unresolved snippets. Existing CI will surface these in build logs.

- **Math formula readability**: Converting `\f$...\f$` to `$...$` assumes the consuming LLM understands LaTeX math notation → Mitigation: all major LLMs handle LaTeX math well. The alternative (stripping math entirely) would lose critical information about what each function computes.
