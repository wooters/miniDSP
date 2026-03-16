## Why

The pre-push hook's version/tag reminder is too terse for AI coding agents. When an agent runs `git push`, it sees a vague "Consider: git tag -a v0.3.0" message but lacks the context to decide whether a version bump is actually needed — it doesn't know what changed since the last tag, whether the VERSION file already matches the tag, or what semantic versioning action to take. This leads to repeated confusion and wasted back-and-forth.

## What Changes

- Enhance the pre-push hook's tag-status section to emit structured, actionable information:
  - Show the latest tag, whether HEAD is already tagged, and the VERSION file value.
  - List a short summary of commits since the last tag (subjects only) so the agent can assess scope.
  - Detect common states (tag matches VERSION, VERSION already bumped but untagged, VERSION not bumped) and print a specific recommended action for each.
  - Use clear labels and machine-friendly formatting so AI agents can parse the status reliably.

## Capabilities

### New Capabilities
- `pre-push-tag-diagnostics`: Structured tag/version diagnostics in the pre-push hook that give AI coding agents enough context to decide whether to bump VERSION and/or create a tag.

### Modified Capabilities

(none — no existing spec-level requirements are changing)

## Impact

- **Code**: `.git/hooks/pre-push` (and the install source at `hooks/pre-push` if one exists)
- **Dependencies**: None — uses only git and shell builtins
- **APIs**: None
- **Systems**: Only affects the local developer workflow; no CI or remote changes
