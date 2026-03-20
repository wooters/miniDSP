## ADDED Requirements

### Requirement: llms.txt index file generation
The build system SHALL generate an `llms.txt` file in `docs/html/` containing a concise markdown index of the miniDSP library. The file SHALL include the library name, brief description, a link to `llms-full.txt`, and a link to the GitHub repository and docs site. The file SHALL be under 1KB.

#### Scenario: llms.txt is generated during docs build
- **WHEN** `make docs` is run
- **THEN** `docs/html/llms.txt` EXISTS and contains the project name "miniDSP", a link to `llms-full.txt`, and a link to the GitHub repo

#### Scenario: llms.txt is deployed to GitHub Pages
- **WHEN** the docs CI workflow completes
- **THEN** `llms.txt` is accessible at the root of the GitHub Pages site (e.g., `https://wooters.github.io/miniDSP/llms.txt`)

---

### Requirement: llms-full.txt complete API reference generation
The build system SHALL generate an `llms-full.txt` file in `docs/html/` containing the complete API reference extracted from Doxygen XML output. For each public function, the file SHALL include the function signature, brief description, detailed description (with math formulas), all `@param` entries, `@return` documentation, and `@code` examples.

#### Scenario: All public functions are documented
- **WHEN** `llms-full.txt` is generated
- **THEN** every non-static function declared in `include/minidsp.h` appears in the output with its signature, parameter docs, and return value description

#### Scenario: Math formulas are converted to standard LaTeX
- **WHEN** Doxygen XML contains formula references (e.g., `<formula>` elements)
- **THEN** the output contains standard LaTeX math delimiters (`$...$` for inline, `$$...$$` for display) instead of Doxygen `\f$...\f$` syntax

#### Scenario: Code examples are included
- **WHEN** a function's Doxygen doc-comment contains an `@code` block
- **THEN** the corresponding code appears in the output as a fenced code block (````c`)

---

### Requirement: Guide content inclusion in llms-full.txt
The build system SHALL append tutorial/guide content from `guides/*.md` files to `llms-full.txt`, with Doxygen-specific syntax stripped and replaced with standard markdown equivalents.

#### Scenario: Doxygen math syntax is converted
- **WHEN** a guide file contains `\f$...\f$` (inline) or `\f[...\f]` (display) math
- **THEN** the output uses `$...$` and `$$...$$` respectively

#### Scenario: HTML-only blocks are removed
- **WHEN** a guide file contains `\htmlonly`...`\endhtmlonly` blocks
- **THEN** those blocks are omitted from the output entirely

#### Scenario: Snippet directives are resolved
- **WHEN** a guide file contains `\snippet filename.c snippet-id`
- **THEN** the output contains the actual code extracted from between `//! [snippet-id]` markers in the referenced source file, rendered as a fenced code block

#### Scenario: Heading anchors are stripped
- **WHEN** a guide file heading contains `{#anchor-id}` suffix
- **THEN** the output contains the heading text without the anchor syntax

#### Scenario: Guide ordering follows tutorials.md
- **WHEN** `llms-full.txt` is generated
- **THEN** guides appear in the order they are listed via `\subpage` directives in `guides/tutorials.md`

---

### Requirement: Doxygen XML output enabled
The Doxyfile SHALL have `GENERATE_XML = YES` so that the generation script has structured API data to parse.

#### Scenario: XML output directory exists after Doxygen runs
- **WHEN** `doxygen Doxyfile` completes
- **THEN** the directory `docs/xml/` EXISTS and contains XML files including `minidsp_8h.xml`

---

### Requirement: Build pipeline integration
The `docs` target in the root Makefile SHALL run the generation script after Doxygen, so that `llms.txt` and `llms-full.txt` are produced as part of every `make docs` invocation. The script SHALL use only Python 3 standard library modules (no pip dependencies).

#### Scenario: make docs produces llms files
- **WHEN** `make docs` is run from the repo root
- **THEN** both `docs/html/llms.txt` and `docs/html/llms-full.txt` exist

#### Scenario: No external Python dependencies
- **WHEN** the generation script is run
- **THEN** it succeeds with only Python 3 standard library (no `pip install` required)

#### Scenario: Script depends on Doxygen XML
- **WHEN** the script is run before `doxygen Doxyfile`
- **THEN** the script exits with a non-zero status and prints an error indicating the XML directory is missing

---

### Requirement: HTML discoverability via link tag
The custom Doxygen header (`doxygen-custom/header.html`) SHALL include a `<link>` tag enabling programmatic discovery of the llms-full.txt file.

#### Scenario: Link tag present in generated HTML
- **WHEN** Doxygen generates HTML output using the custom header
- **THEN** every HTML page contains `<link rel="help" type="text/markdown" href="llms-full.txt">` in the `<head>` section

---

### Requirement: Generated files are not tracked in git
The `llms.txt` and `llms-full.txt` files in `docs/html/` SHALL NOT be tracked in version control. They are build artifacts generated on every docs build.

#### Scenario: Files excluded from git
- **WHEN** `make docs` is run and `docs/html/llms.txt` and `docs/html/llms-full.txt` are generated
- **THEN** `git status` does not show them as untracked (they are covered by existing `docs/html/` gitignore or an explicit entry)

---

### Requirement: Clean target removes XML output
The `make clean` target SHALL remove the `docs/xml/` directory along with `docs/html/`.

#### Scenario: Clean removes XML
- **WHEN** `make clean` is run
- **THEN** `docs/xml/` does not exist
