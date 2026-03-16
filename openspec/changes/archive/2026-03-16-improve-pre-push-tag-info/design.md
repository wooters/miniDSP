## Context

The pre-push hook (`scripts/pre-push`, installed to `.git/hooks/pre-push`) runs tests before allowing pushes to main, then prints a version/tag reminder. The current reminder is a 4-line block that says how many commits exist since the last tag and suggests a `git tag` command. AI coding agents lack the context to act on this: they don't know whether the VERSION file was already bumped, whether commits are trivial or substantial, or what the correct next step is.

The hook source lives at `scripts/pre-push` and is copied to `.git/hooks/pre-push` via `make install-hooks`.

## Goals / Non-Goals

**Goals:**
- Emit enough structured information that an AI agent can determine the correct action without further investigation
- Detect the three common states: (1) HEAD already tagged, (2) VERSION bumped but tag not yet created, (3) VERSION not bumped since last tag
- Print a specific, copy-pasteable recommended action for each state
- Include a commit summary so the agent can assess whether changes warrant a patch, minor, or major bump

**Non-Goals:**
- Automatically bumping VERSION or creating tags (the hook only advises)
- Changing the test-gating behavior of the hook
- Parsing commit messages for conventional-commit semantics
- Supporting non-main branches

## Decisions

**Decision 1: Replace the existing version-check block with a multi-state diagnostic**

The current code has a single code path (HEAD untagged → print reminder). The new code will detect three states and print tailored guidance for each:

1. **HEAD is already tagged** → print "TAG OK" confirmation with the tag name. No action needed.
2. **VERSION > latest tag** (bumped but untagged) → print "READY TO TAG" with the exact `git tag` command.
3. **VERSION == latest tag version** (not bumped) → print "VERSION NOT BUMPED" with commit summary and guidance to update VERSION first, then tag.

*Alternative considered*: Outputting machine-parseable JSON. Rejected because the hook output appears in terminal alongside other git output, and shell-friendly labeled lines are easier for both humans and AI agents to read.

**Decision 2: Show up to 10 commit subjects since last tag**

When VERSION hasn't been bumped, the agent needs to see what changed to decide on patch vs minor vs major. `git log --oneline` with a limit of 10 provides enough context without flooding the terminal. If there are more than 10, show the count.

*Alternative considered*: Showing full commit messages. Rejected — subjects are sufficient for triage and keep output compact.

**Decision 3: Compare version semantically using string equality**

Compare the VERSION file content (e.g., `0.3.0`) against the latest tag stripped of its `v` prefix. Simple string equality works because the project uses strict `vMAJOR.MINOR.PATCH` tags and bare `MAJOR.MINOR.PATCH` in VERSION.

## Risks / Trade-offs

- **[Risk] Hook output becomes verbose** → Mitigated by using clear section headers and keeping commit list capped at 10 lines. The "TAG OK" path is actually shorter than the current output.
- **[Risk] Edge case: no tags exist at all** → Keep existing behavior (skip the diagnostic entirely when `git describe` finds nothing).
- **[Risk] VERSION file missing or empty** → Already handled by the existing `${version:-unknown}` fallback; preserve this.
