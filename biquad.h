#ifndef BIQUAD_H
#define BIQUAD_H

#include <math.h>
#include <stdlib.h>

#ifndef M_LN2
#define M_LN2 0.69314718055994530942
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* whatever sample type you want */
typedef double smp_type;

/* this holds the data required to update samples thru a filter */
typedef struct {
    smp_type a0, a1, a2, a3, a4;
    smp_type x1, x2, y1, y2;
}
biquad;

extern smp_type BiQuad(const smp_type sample, biquad* const b);
extern biquad *BiQuad_new(const int type, smp_type dbGain, /* gain of filter */
                          const smp_type freq,             /* center frequency */
                          const smp_type srate,            /* sampling rate */
                          const smp_type bandwidth);       /* bandwidth in octaves */

/* filter types */
enum FILT_TYPE {
    LPF, /* low pass filter */
    HPF, /* High pass filter */
    BPF, /* band pass filter */
    NOTCH, /* Notch Filter */
    PEQ, /* Peaking band EQ filter */
    LSH, /* Low shelf filter */
    HSH /* High shelf filter */
};

#endif
