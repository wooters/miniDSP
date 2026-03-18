## ADDED Requirements

### Requirement: Guide page exists at guides/resample-tool.md
A Doxygen guide page SHALL exist at `guides/resample-tool.md` with the anchor `{#resample-tool}`.

#### Scenario: Page renders in Doxygen
- **WHEN** `make docs` is run
- **THEN** the resample tool guide appears in the generated HTML documentation

### Requirement: Guide documents CLI usage
The guide SHALL document the tool's command-line interface: positional arguments (input file, target rate, output file) and optional flags (`-z`, `-b`).

#### Scenario: Usage section present
- **WHEN** a user reads the guide
- **THEN** they see the full usage syntax, a description of each argument, and at least one example command

### Requirement: Guide explains -z and -b options
The guide SHALL explain what the `-z` (zero-crossings) and `-b` (kaiser beta) options control, including their defaults and effect on output quality.

#### Scenario: Options documented with defaults
- **WHEN** a user reads the options section
- **THEN** they see that `-z` defaults to 32 and `-b` defaults to 10.0, with an explanation of what each controls

### Requirement: Guide links to resampling API page
The guide SHALL link to the existing `resampling` guide page for users who want deeper mathematical detail on the polyphase sinc resampler.

#### Scenario: Cross-reference present
- **WHEN** a user wants to understand the math
- **THEN** the guide provides a link to `\ref resampling`

### Requirement: Guide is linked from tutorials.md
The guide SHALL be reachable via a `\subpage resample-tool` entry in `guides/tutorials.md` under the Tools section.

#### Scenario: Navigation link present
- **WHEN** a user browses the Tutorials page
- **THEN** they see a link to the resample tool guide under Tools

### Requirement: Doxyfile EXAMPLE_PATH includes tools/resample
The Doxyfile `EXAMPLE_PATH` SHALL include `tools/resample` so that `\snippet resample.c` directives work.

#### Scenario: Snippet directive resolves
- **WHEN** the guide uses `\snippet resample.c` with a valid marker
- **THEN** Doxygen resolves it and embeds the code
