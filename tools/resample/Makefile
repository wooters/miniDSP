
# Makefile for resample — sample rate conversion tool
#
# Build:    make
# Clean:    make clean

include ../../config.mk

CFLAGS += -O2 -g -I../../include

default: resample

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../../libminidsp.a
../../libminidsp.a:
	$(MAKE) -C ../../ libminidsp.a

resample: resample.c ../../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L../../ -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

.PHONY: clean
clean:
	-rm -rf resample resample.dSYM
