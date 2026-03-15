/**
 * @file test_helpers.h
 * @brief Shared test infrastructure for the miniDSP test suite.
 *
 * Provides global test counters, the RUN_TEST macro, and common helpers
 * used across all per-module test files.
 */

#ifndef TEST_HELPERS_H
#define TEST_HELPERS_H

#include <stdio.h>
#include <math.h>

/* -----------------------------------------------------------------------
 * Global test counters — defined in test_minidsp.c (the driver).
 * -----------------------------------------------------------------------*/

extern int tests_run;
extern int tests_passed;
extern int tests_failed;

/* -----------------------------------------------------------------------
 * Helpers
 * -----------------------------------------------------------------------*/

/** Compare two doubles with a tolerance.  Returns 1 if close enough. */
static inline int approx_equal(double a, double b, double tolerance)
{
    return fabs(a - b) <= tolerance;
}

/** Run a single test and print the result. */
#define RUN_TEST(test_func)                                       \
    do {                                                          \
        tests_run++;                                              \
        printf("  %-50s", #test_func);                            \
        if (test_func()) {                                        \
            tests_passed++;                                       \
            printf("[PASS]\n");                                   \
        } else {                                                  \
            tests_failed++;                                       \
            printf("[FAIL]\n");                                   \
        }                                                         \
    } while (0)

/* -----------------------------------------------------------------------
 * Helper: create a delayed copy of a signal using circular rotation.
 *
 * This mimics what happens when a sound reaches one microphone later
 * than another.  A positive 'delay' shifts the signal to the right
 * (later in time).
 * -----------------------------------------------------------------------*/
static inline void delay_signal(const double *in, double *out,
                                unsigned n, int delay)
{
    for (unsigned i = 0; i < n; i++) {
        /* Use modular arithmetic with unsigned wraparound.
         * For negative delays, the unsigned conversion + modulus
         * produces the correct circular shift. */
        unsigned j = (i + (unsigned)delay) % n;
        out[j] = in[i];
    }
}

#endif /* TEST_HELPERS_H */
