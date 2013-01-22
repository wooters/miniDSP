#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <inttypes.h>
#include <sndfile.h>


typedef struct HTKheader {
  uint32_t nvecs;
  uint32_t sampperiod;
  uint16_t vecsize;
  uint16_t parmkind;
} HTKheader;

void SwapBytes(void* pv, size_t n);

/* Read an audio file using libsndfile */
void  read_audio(const char* infile, 
		 float** indata, 
		 size_t* datalen, 
		 unsigned* samprate, 
		 unsigned donorm);

/* Write the feature vectors in HTK feature file format */
void write_feats(const char* outfile, 
		 float** outvecs, 
		 size_t nvecs, 
		 size_t veclen, 
		 unsigned vecsamprate);

/* Convenience routine to allocate memory */
void* malloc_or_die(size_t nbytes, const char* msg);

#endif
