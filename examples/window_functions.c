/**
 * @file window_functions.c
 * @brief Example: generate and compare common FFT window functions.
 *
 * Build and run (from repository root):
 *   make
 *   make -C examples window_functions
 *   ./examples/window_functions
 */

#include <stdio.h>
#include "minidsp.h"

static void print_summary(const char *name, const double *w, unsigned n)
{
    printf("%-10s first=% .6f  center=% .6f  last=% .6f\n",
           name, w[0], w[n / 2], w[n - 1]);
}

int main(void)
{
    const unsigned N = 64;
    double hann[N], hamming[N], blackman[N], rect[N];

    //! [generate-hann]
    MD_Gen_Hann_Win(hann, N);
    //! [generate-hann]

    //! [generate-hamming]
    MD_Gen_Hamming_Win(hamming, N);
    //! [generate-hamming]

    //! [generate-blackman]
    MD_Gen_Blackman_Win(blackman, N);
    //! [generate-blackman]

    //! [generate-rect]
    MD_Gen_Rect_Win(rect, N);
    //! [generate-rect]

    //! [generate-all-windows]
    MD_Gen_Hann_Win(hann, N);
    MD_Gen_Hamming_Win(hamming, N);
    MD_Gen_Blackman_Win(blackman, N);
    MD_Gen_Rect_Win(rect, N);
    //! [generate-all-windows]

    printf("Window summaries (N=%u)\n", N);
    print_summary("Hanning", hann, N);
    print_summary("Hamming", hamming, N);
    print_summary("Blackman", blackman, N);
    print_summary("Rect", rect, N);

    return 0;
}
