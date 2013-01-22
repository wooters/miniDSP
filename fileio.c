
#include "utils.h"

const int _i = 1;
#define is_bigendian() ( (*(char*)&_i) == 0 )

void SwapBytes(void* pv, size_t n) {
  char *p = pv;
  size_t lo, hi;
  for(lo=0, hi=n-1; hi>lo; lo++, hi--) {
    char tmp=p[lo];
    p[lo] = p[hi];
    p[hi] = tmp;
  }
}
#define SWAP(x) SwapBytes(&x, sizeof(x));

/* Function to read an audio file */
void  read_audio(const char* infile, 
		 float** indata, 
		 size_t* datalen, 
		 unsigned* samprate, 
		 unsigned donorm) {

  /* The libsndfile file structure */
  SF_INFO sfinfo;

  /* Open and read the audio file info */
  SNDFILE* sf = sf_open(infile, SFM_READ, &sfinfo);
  if (sfinfo.channels != 1) {
    printf("Input file has %u channels. We can only handle 1-channel data\n",
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
  float* tmpdata = malloc_or_die(nsamps*sizeof(float), 
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

void write_feats(const char* outfile, 
		 float** outvecs, 
		 size_t nvecs, 
		 size_t veclen, 
		 unsigned vecsamprate) {

  HTKheader hdr;
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
  fwrite(&hdr,sizeof(HTKheader),1,f);

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

/* these two pragma's stop the warnings that gcc gives */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
#pragma GCC diagnostic ignored "-Wformat-security"
void* malloc_or_die(size_t nbytes, const char* msg) {

    void* tmp = malloc(nbytes);
    if (tmp == NULL) {
      printf(msg);
      exit(1);
    }
    return tmp;
}
#pragma GCC diagnostic pop



