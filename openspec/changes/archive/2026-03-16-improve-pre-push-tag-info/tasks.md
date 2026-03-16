## 1. Implement tag diagnostic logic

- [x] 1.1 Replace the existing version-check block in `scripts/pre-push` (lines 46-60) with the new multi-state diagnostic that detects three states: HEAD tagged, VERSION bumped but untagged, VERSION not bumped
- [x] 1.2 In the HEAD-already-tagged path, print "TAG OK" with the tag name(s) and skip further diagnostics
- [x] 1.3 In the VERSION-bumped-but-untagged path, print "READY TO TAG" with the exact `git tag -a` command
- [x] 1.4 In the VERSION-not-bumped path, print "VERSION NOT BUMPED" with commit count, up to 10 commit subjects, and guidance to update VERSION first
- [x] 1.5 Always show labeled state values (latest tag, VERSION file, commits since tag) in the diagnostic section

## 2. Install and verify

- [x] 2.1 Run `make install-hooks` to copy updated script to `.git/hooks/pre-push`
- [x] 2.2 Manually test each state: verify output when HEAD is tagged, when VERSION is bumped but untagged, and when VERSION matches the latest tag with new commits
