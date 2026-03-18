# Error Handling {#error-handling}

miniDSP uses a configurable error handler to report precondition
violations.  The library **never** calls `abort()` or terminates the
host process.  When a function detects bad input (NULL pointer, invalid
size, out-of-range parameter), it:

1. Calls the active error handler with an error code, function name,
   and human-readable message.
2. Returns a **safe default** value appropriate to the function's
   return type.

## Error codes

| Code | Meaning |
|------|---------|
| `MD_ERR_NULL_POINTER`  | A required pointer argument is NULL |
| `MD_ERR_INVALID_SIZE`  | A size or count argument is invalid (e.g. N == 0) |
| `MD_ERR_INVALID_RANGE` | A range or bound is invalid (e.g. min >= max) |
| `MD_ERR_ALLOC_FAILED`  | A memory allocation failed |

## Default behavior

If no custom handler is installed, errors are printed to `stderr`:

```
miniDSP: error 1 in MD_energy(): a is NULL
```

The function then returns a safe default:
- `void` functions: early return (no-op)
- `double`-returning functions: `0.0`
- `unsigned`-returning functions: `0`
- `int`-returning functions: `-1`

## Installing a custom handler

Use MD_set_error_handler() to install your own handler.  The handler
receives the error code, function name, and a descriptive message.

```c
static void my_handler(MD_ErrorCode code, const char *fn, const char *msg) {
    app_log(LOG_WARN, "miniDSP error %d in %s: %s", code, fn, msg);
}

int main(void) {
    MD_set_error_handler(my_handler);
    /* ... use miniDSP ... */
}
```

Pass `NULL` to restore the default stderr handler:

```c
MD_set_error_handler(NULL);  /* back to default */
```

## Thread safety

MD_set_error_handler() must be called **once, before any other
`MD_*` function**.  It must not be called concurrently with any other
miniDSP function.  This is the same threading contract as FFTW's
initialization functions.

## Sentinel returns vs errors

Some functions return sentinel values for valid-but-unresolved runtime
outcomes.  For example, MD_f0_autocorrelation() returns `0.0` when no
pitch peak is found.  These are **not errors** and do not trigger the
error handler.  The error handler fires only for **precondition
violations** (bad arguments that represent programmer mistakes).
