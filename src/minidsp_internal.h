/**
 * @file minidsp_internal.h
 * @brief Internal header for cross-file dependencies within the minidsp module.
 *
 * This header is NOT part of the public API.  It declares symbols shared
 * between the split minidsp_*.c translation units but not exposed to
 * library consumers.
 */

#ifndef MINIDSP_INTERNAL_H
#define MINIDSP_INTERNAL_H

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Error reporting (see minidsp_error.c)
 * -----------------------------------------------------------------------*/

/** Report a precondition violation to the active error handler. */
void md_report_error(MD_ErrorCode code, const char *func_name,
                     const char *message);

/**
 * Check a precondition in a function that returns a value.
 * On failure: report the error, then return @a retval.
 */
#define MD_CHECK(cond, code, msg, retval) do { \
    if (!(cond)) {                              \
        md_report_error((code), __func__, (msg)); \
        return (retval);                        \
    }                                           \
} while (0)

/**
 * Check a precondition in a void function.
 * On failure: report the error, then return.
 */
#define MD_CHECK_VOID(cond, code, msg) do { \
    if (!(cond)) {                           \
        md_report_error((code), __func__, (msg)); \
        return;                              \
    }                                        \
} while (0)

/* -----------------------------------------------------------------------
 * Cross-file teardown hooks (called from MD_shutdown)
 * -----------------------------------------------------------------------*/

/** Tear down the GCC-PHAT FFT cache and reset its cached length.
 *  Called only from MD_shutdown(). */
void md_gcc_teardown(void);

/** Tear down the BFSK sine carrier cache.
 *  Called only from MD_shutdown(). */
void md_steg_teardown(void);

#endif /* MINIDSP_INTERNAL_H */
