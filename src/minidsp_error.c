/**
 * @file minidsp_error.c
 * @brief Configurable error handler for precondition violations.
 * @author Chuck Wooters <wooters@hey.com>
 */

#include "minidsp.h"

/* -----------------------------------------------------------------------
 * Default handler: log to stderr, never abort.
 * -----------------------------------------------------------------------*/

static void default_error_handler(MD_ErrorCode code,
                                  const char  *func_name,
                                  const char  *message)
{
    fprintf(stderr, "miniDSP: error %d in %s(): %s\n",
            (int)code, func_name, message);
}

/* -----------------------------------------------------------------------
 * Handler state — set once before use, not thread-safe.
 * -----------------------------------------------------------------------*/

static MD_ErrorHandler current_handler = default_error_handler;

void MD_set_error_handler(MD_ErrorHandler handler)
{
    current_handler = handler ? handler : default_error_handler;
}

/* -----------------------------------------------------------------------
 * Internal reporting function (called by MD_CHECK macros).
 * -----------------------------------------------------------------------*/

void md_report_error(MD_ErrorCode code, const char *func_name,
                     const char *message)
{
    current_handler(code, func_name, message);
}
