
# Makefile for audio_steg — audio steganography tool
#
# Build:    make
# Clean:    make clean

include ../../config.mk

CFLAGS += -O2 -g -I../../include

default: audio_steg

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../../libminidsp.a
../../libminidsp.a:
	$(MAKE) -C ../../ libminidsp.a

audio_steg: audio_steg.c ../../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L../../ -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

.PHONY: clean
clean:
	-rm -rf audio_steg audio_steg.dSYM
