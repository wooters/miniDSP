/**
 * @file fileio.h
 *
 */
#ifndef FILEIO_H
#define FILEIO_H

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <inttypes.h>
#include <sndfile.h>


/* Read an audio file using libsndfile */
void  FIO_read_audio(const char* infile, 
		     float** indata, 
		     size_t* datalen, 
		     unsigned* samprate, 
		     unsigned donorm);

/* Write feature vectors in HTK feature file format */
void FIO_write_htk_feats(const char* outfile, 
			 const float** const outvecs, 
			 const size_t nvecs, 
			 const size_t veclen, 
			 const unsigned vecsamprate);


#endif
