/**
 * @file test_minidsp.c
 * @brief Test driver for the miniDSP library.
 *
 * This slim driver defines the global test counters, calls each module's
 * run_*_tests() function, prints the summary, and calls MD_shutdown().
 *
 * How to compile (from the tests/ directory):
 *   make test_minidsp
 *
 * How to run:
 *   ./test_minidsp
 */

#include <stdio.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * Global test counters (declared extern in test_helpers.h)
 * -----------------------------------------------------------------------*/

int tests_run    = 0;
int tests_passed = 0;
int tests_failed = 0;

/* -----------------------------------------------------------------------
 * Forward declarations for per-module test runners
 * -----------------------------------------------------------------------*/

void run_core_tests(void);
void run_effects_tests(void);
void run_generators_tests(void);
void run_spectrum_tests(void);
void run_fir_tests(void);
void run_gcc_tests(void);
void run_biquad_tests(void);
void run_windows_tests(void);
void run_dtmf_tests(void);
void run_spectext_tests(void);
void run_steg_tests(void);
void run_fileio_tests(void);
void run_resample_tests(void);

/* -----------------------------------------------------------------------
 * Main: run all tests
 * -----------------------------------------------------------------------*/

int main(void)
{
    printf("=== miniDSP Test Suite ===\n\n");

    run_core_tests();
    run_effects_tests();
    run_fir_tests();
    run_windows_tests();
    run_spectrum_tests();
    run_gcc_tests();
    run_biquad_tests();
    run_generators_tests();
    run_dtmf_tests();
    run_spectext_tests();
    run_steg_tests();
    run_fileio_tests();
    run_resample_tests();

    /* Clean up FFTW resources */
    MD_shutdown();

    /* Print summary */
    printf("\n=== Results: %d/%d passed", tests_passed, tests_run);
    if (tests_failed > 0) {
        printf(", %d FAILED", tests_failed);
    }
    printf(" ===\n");

    return (tests_failed > 0) ? 1 : 0;
}
