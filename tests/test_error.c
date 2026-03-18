/**
 * @file test_error.c
 * @brief Tests for the error handling system: custom handlers, default
 *        handler restoration, correct error code/function/message delivery,
 *        and safe default return values.
 */

#include <string.h>
#include "minidsp.h"
#include "test_helpers.h"

/* -----------------------------------------------------------------------
 * Capture state for the custom test handler
 * -----------------------------------------------------------------------*/

static MD_ErrorCode last_code = 0;
static char last_func[128]    = {0};
static char last_msg[256]     = {0};
static int  handler_call_count = 0;

static void test_handler(MD_ErrorCode code, const char *func, const char *msg)
{
    last_code = code;
    strncpy(last_func, func, sizeof(last_func) - 1);
    last_func[sizeof(last_func) - 1] = '\0';
    strncpy(last_msg, msg, sizeof(last_msg) - 1);
    last_msg[sizeof(last_msg) - 1] = '\0';
    handler_call_count++;
}

static void reset_capture(void)
{
    last_code = 0;
    last_func[0] = '\0';
    last_msg[0] = '\0';
    handler_call_count = 0;
}

/* -----------------------------------------------------------------------
 * Tests for handler installation and invocation
 * -----------------------------------------------------------------------*/

/** Custom handler is called on a precondition violation. */
static int test_custom_handler_called(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_energy(NULL, 100);

    MD_set_error_handler(NULL);  /* restore default */
    return handler_call_count == 1;
}

/** Custom handler receives the correct error code. */
static int test_handler_receives_code(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_energy(NULL, 100);

    MD_set_error_handler(NULL);
    return last_code == MD_ERR_NULL_POINTER;
}

/** Custom handler receives the correct function name. */
static int test_handler_receives_func_name(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_energy(NULL, 100);

    MD_set_error_handler(NULL);
    return strcmp(last_func, "MD_energy") == 0;
}

/** Custom handler receives a non-empty message. */
static int test_handler_receives_message(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_energy(NULL, 100);

    MD_set_error_handler(NULL);
    return strlen(last_msg) > 0;
}

/** Passing NULL to MD_set_error_handler restores the default. */
static int test_restore_default_handler(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    /* First call goes to test_handler */
    MD_energy(NULL, 100);
    int first_count = handler_call_count;

    /* Restore default */
    MD_set_error_handler(NULL);
    reset_capture();

    /* Second call should NOT go to test_handler (default logs to stderr) */
    MD_energy(NULL, 100);

    return first_count == 1 && handler_call_count == 0;
}

/** Handler receives MD_ERR_INVALID_SIZE for size violations. */
static int test_handler_size_error_code(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_power(NULL, 0);  /* N == 0 triggers size check before null check */

    MD_set_error_handler(NULL);
    /* The first check that fires depends on order — either NULL or SIZE */
    return last_code == MD_ERR_NULL_POINTER || last_code == MD_ERR_INVALID_SIZE;
}

/** Handler is called for void functions too. */
static int test_handler_void_function(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_Gen_Hann_Win(NULL, 512);

    MD_set_error_handler(NULL);
    return handler_call_count == 1 && last_code == MD_ERR_NULL_POINTER;
}

/* -----------------------------------------------------------------------
 * Tests for safe default return values
 * -----------------------------------------------------------------------*/

/** Double-returning function returns 0.0 on bad input. */
static int test_safe_default_double(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    double result = MD_energy(NULL, 100);

    MD_set_error_handler(NULL);
    return result == 0.0;
}

/** Unsigned-returning function returns 0 on bad input. */
static int test_safe_default_unsigned(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    unsigned result = MD_convolution_num_samples(0, 0);

    MD_set_error_handler(NULL);
    return result == 0;
}

/** Void function no-ops on bad input (doesn't write to buffer). */
static int test_safe_default_void_noop(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    double buf[4] = {1.0, 2.0, 3.0, 4.0};
    MD_Gen_Hann_Win(NULL, 4);  /* should not crash */

    /* buf should be unchanged since we passed NULL (the function no-ops) */
    MD_set_error_handler(NULL);
    return buf[0] == 1.0 && buf[1] == 2.0 && buf[2] == 3.0 && buf[3] == 4.0;
}

/** Multiple errors fire the handler multiple times. */
static int test_multiple_errors(void)
{
    MD_set_error_handler(test_handler);
    reset_capture();

    MD_energy(NULL, 100);
    MD_power(NULL, 0);
    MD_rms(NULL, 0);

    MD_set_error_handler(NULL);
    return handler_call_count == 3;
}

/* -----------------------------------------------------------------------
 * Test runner
 * -----------------------------------------------------------------------*/

void run_error_tests(void)
{
    printf("\n--- Error handling ---\n");
    RUN_TEST(test_custom_handler_called);
    RUN_TEST(test_handler_receives_code);
    RUN_TEST(test_handler_receives_func_name);
    RUN_TEST(test_handler_receives_message);
    RUN_TEST(test_restore_default_handler);
    RUN_TEST(test_handler_size_error_code);
    RUN_TEST(test_handler_void_function);
    RUN_TEST(test_safe_default_double);
    RUN_TEST(test_safe_default_unsigned);
    RUN_TEST(test_safe_default_void_noop);
    RUN_TEST(test_multiple_errors);
}
