## Context

The README's audio I/O bullet currently reads:

> Read audio files in any format supported by libsndfile (WAV, FLAC, AIFF, OGG, etc.)

libsndfile actually supports 20+ major container formats. Adding a count and a documentation link helps users understand the breadth without needing to look it up themselves.

## Goals / Non-Goals

**Goals:**
- Add an approximate format count (~20+) to the existing README bullet.
- Add a hyperlink to the official libsndfile formats page.

**Non-Goals:**
- Listing every individual format in the README.
- Changing any code or library behavior.

## Decisions

- **Say "20+" rather than an exact count.** libsndfile's format list grows across versions (MP3 was added in 1.1.0). An approximate number stays accurate longer.
- **Link to `https://libsndfile.github.io/libsndfile/formats.html`** — the canonical formats page hosted by the libsndfile project.

## Risks / Trade-offs

- The linked page could move if libsndfile reorganizes their docs. Low risk — the URL has been stable for years.
