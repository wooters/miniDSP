
# Makefile for mel_viz — mel-spectrum audio visualizer
#
# Build:    make
# Clean:    make clean

include ../../config.mk

CFLAGS += -O2 -g -I../../include

default: mel_viz

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../../libminidsp.a
../../libminidsp.a:
	$(MAKE) -C ../../ libminidsp.a

mel_viz: mel_viz.c ../../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L../../ -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

.PHONY: clean
clean:
	-rm -rf mel_viz mel_viz.dSYM
