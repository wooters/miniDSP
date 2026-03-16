## 1. Platform compatibility note

- [x] 1.1 Add a sentence after line 8 ("A small C library...") stating the library compiles and runs on Ubuntu and macOS, that Windows is untested, and PRs are welcome

## 2. Section reorder

- [x] 2.1 Move the entire "Build and Test" section (lines 423–474: heading, Dependencies, Compile, Test, Container, Docs, Hooks) to appear before "Use in your project"

## 3. Fix Apple container dependency

- [x] 3.1 In the Dependencies table, change the Apple container macOS column from "macOS 26+ built-in" to "Install from [GitHub](https://github.com/apple/container)"

## 4. Add audio_steg tool

- [x] 4.1 Add an `audio_steg` entry to the "Tools" section (after `mel_viz`) with a heading, one-line description, and build/run snippet matching the `mel_viz` style
