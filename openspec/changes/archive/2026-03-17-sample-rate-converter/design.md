## Context

miniDSP provides `MD_resample()` (polyphase sinc interpolation) and file I/O via `FIO_read_audio()` / `FIO_write_wav()`. No CLI tool wraps these together. The existing tools (`mel_viz`, `audio_steg`) in `tools/` demonstrate the pattern: a standalone C program with its own `Makefile` that links `libminidsp.a`.

Key constraint: `FIO_read_audio()` returns `float*` but `MD_resample()` operates on `double*`, so the tool must convert between the two.

## Goals / Non-Goals

**Goals:**
- Provide a single-command sample rate conversion: `resample <input> <rate> <output>`
- Use sensible resampler defaults (32 zero-crossings, kaiser beta 10.0) while allowing optional overrides
- Print input/output metadata (sample rate, duration, sample count) for user feedback
- Validate constraints up front with clear error messages

**Non-Goals:**
- Multi-channel support (the library only handles mono)
- Output formats other than WAV (libsndfile write path is WAV-only)
- Real-time / streaming resampling
- GUI or interactive mode

## Decisions

**1. CLI argument structure**: positional args for the three required values, optional flags for resampler quality.

```
resample <input-file> <target-rate> <output.wav>
resample -z 64 -b 14.0 <input-file> <target-rate> <output.wav>
```

*Rationale*: Positional args keep the common case simple. Optional `-z` (zero-crossings) and `-b` (kaiser beta) flags let advanced users tune quality without cluttering the default invocation. Alternative considered: all flags (`-i`, `-r`, `-o`) — rejected as unnecessarily verbose for a three-argument tool.

**2. Float/double conversion**: Convert `float*` from `FIO_read_audio()` to `double*` for `MD_resample()`, then convert `double*` output back to `float*` for `FIO_write_wav()`.

*Rationale*: The library's DSP pipeline works in double precision. The file I/O layer uses float (libsndfile convention). Two simple loops handle the conversion. Alternative considered: modifying `FIO_read_audio` to return doubles — rejected as it would change the public API.

**3. Memory allocation**: `malloc` all buffers (input float, input double, output double, output float). Use `MD_resample_output_len()` to size the output buffer. Free each buffer as soon as it is no longer needed to minimize peak memory:

1. Read input → `float*` allocated by `FIO_read_audio()`
2. Convert to `double*` → free the input `float*`
3. Resample → free the input `double*`
4. Convert output to `float*` → free the output `double*`
5. Write WAV → free the output `float*`

*Rationale*: Audio files can be arbitrarily large; stack allocation is not safe. Eagerly freeing buffers keeps peak memory at roughly 2x the larger of the input or output signal, rather than holding all four buffers simultaneously.

**4. Argument parsing**: Hand-rolled `getopt`-style loop for `-z` and `-b`, then consume remaining positional args.

*Rationale*: Only two optional flags — a full argument parsing library is overkill. `getopt()` is POSIX but the codebase avoids `_XOPEN_SOURCE`. A simple manual loop over `argv` is sufficient and keeps the code self-contained.

## Risks / Trade-offs

- **Large file memory usage** — The tool loads the entire file into memory. For very long recordings this could be significant. → Acceptable for an offline batch tool; streaming would add complexity with no clear user need.
- **No-op resampling** — If input rate equals target rate, the tool still runs the resampler. → Detect this case and copy directly, printing a warning.
- **Non-WAV output extension** — User might specify `output.flac` but get WAV content. → Warn if the output filename doesn't end in `.wav`, but proceed (the data is still valid WAV).
