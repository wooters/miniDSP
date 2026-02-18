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

/** Tear down the GCC-PHAT FFT cache and reset its cached length.
 *  Called only from MD_shutdown(). */
void md_gcc_teardown(void);

#endif /* MINIDSP_INTERNAL_H */
