/**
 * @file fileio.c
 */
#include "fileio.h"

typedef struct FIO_HTKheader {
  uint32_t nvecs;
  uint32_t sampperiod;
  uint16_t vecsize;
  uint16_t parmkind;
} FIO_HTKheader;

const int _i = 1;
#define is_bigendian() ( (*(char*)&_i) == 0 )

void _SwapBytes(void* pv, size_t n) {
  char *p = pv;
  size_t lo, hi;
  for(lo=0, hi=n-1; hi>lo; lo++, hi--) {
    char tmp=p[lo];
    p[lo] = p[hi];
    p[hi] = tmp;
  }
}
#define SWAP(x) _SwapBytes(&x, sizeof(x));

/* these two pragma's stop the warnings that gcc gives */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#pragma GCC diagnostic ignored "-Wformat-security"
void* _malloc_or_die(size_t nbytes, const char* const msg) {

    void* tmp = malloc(nbytes);
    if (tmp == NULL) {
      printf(msg);
      exit(1);
    }
    return tmp;
}
#pragma GCC diagnostic pop


/** 
 * Read an audio file from disk.
 *
 * @param infile Then name of the file to read from.
 *
 * @param indata A pointer to a float pointer. The address of the
 * read-in audio data will be stored at this location.
 *
 * @param datalen A pointer to a size_t variable in which to store the
 * length (in samples) of the input data.
 *
 * @param samprate A pointer to an unsigned variable in which to store
 * the sampling rate of the audio data.
 *
 * @param donorm Either a 1 or 0 to indicate if the audio should be
 * normalized when read. If \a donorm == 1, then the audio data will
 * have a range of [-1.0, 1.0]
 *
 * @note Currently only handles single-channel audio files. You can use the 'sox' 
 * utility to split a multi-channel file into multiple single-channel files.
 */
void  FIO_read_audio(const char* const infile, 
		     float** indata, 
		     size_t* datalen, 
		     unsigned* samprate, 
		     unsigned donorm) {

  /* The libsndfile file structure */
  SF_INFO sfinfo;

  /* Open and read the audio file info */
  SNDFILE* sf = sf_open(infile, SFM_READ, &sfinfo);
  if (sfinfo.channels != 1) {
    printf("Input file has %u channels. We only handle 1-channel data at this point\n",
	   sfinfo.channels);
    exit(1);
  }
  sf_count_t nsamps = sfinfo.frames;

  /* Set normalization of sample values on/off */
  if (donorm == 1)
    sf_command(sf, SFC_SET_NORM_FLOAT, NULL, SF_TRUE);
  else
    sf_command(sf, SFC_SET_NORM_FLOAT, NULL, SF_FALSE);

  /* Allocate space for the audio data and read it in */
  float* tmpdata = _malloc_or_die(nsamps*sizeof(float), 
				  "Error allocating memory for audio data\n");
  sf_count_t nread  = sf_read_float(sf, tmpdata, nsamps);
  if (nread != nsamps) {
    printf("Error reading from %s, expected to read %d, instead read %d\n",
	   infile, (int)nsamps, (int)nread);
    exit(1);
  }

  /* Close the sound file */
  sf_close(sf);

  /* Set the output variables */
  *indata   = tmpdata;
  *datalen  = (size_t) nsamps;
  *samprate = (unsigned) sfinfo.samplerate;
}

/**
 * Write an HTK feature file.
 */
void FIO_write_htk_feats(const char* outfile, 
			 const float** const outvecs,
			 const size_t nvecs, 
			 const size_t veclen, 
			 const unsigned vecsamprate) {

  FIO_HTKheader hdr;
  hdr.nvecs = nvecs;
  /* HTK uses 100ns units for the sampling period */
  hdr.sampperiod = (uint32_t)(1.0 / (float)vecsamprate * 1e7);
  hdr.vecsize = veclen*sizeof(float);
  hdr.parmkind = 9; /* 9 == USER (i.e. user-defined type) */

  FILE* f = fopen(outfile, "wb");

  /* 
   * Default byte ordering in HTK files is big-endian, so swap ordering
   * if needed before writing the header.
   */
  if (!is_bigendian()) {  
    SWAP(hdr.nvecs);
    SWAP(hdr.sampperiod);
    SWAP(hdr.vecsize);
    SWAP(hdr.parmkind);
  }
  /* Write the 12-byte HTK header */
  fwrite(&hdr,sizeof(FIO_HTKheader),1,f);

  /* Write the feature vectors */
  for(size_t i=0; i<nvecs; i++) {
    if (!is_bigendian()) {
      for(size_t j=0; j<veclen; j++) {
	float tmp = outvecs[i][j];
	SWAP(tmp);
	fwrite(&tmp,sizeof(float),1,f);
      }
    } else {
      fwrite(outvecs[i],sizeof(float),veclen,f);
    }
  }
  fclose(f);
}



