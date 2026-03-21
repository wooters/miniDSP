
# Test Makefile for miniDSP
#
# Build the comprehensive test suite:
#   make test_minidsp
#
# Run the tests:
#   make test
#
# Build the legacy demo programs (require gnuplot_i / portaudio):
#   make gcc_phat_test
#   make gcc_phat_multi_test

include ../config.mk

CFLAGS += -g -I../include

# Paths to optional dependencies (legacy demos)
GP = $(HOME)/src/gnuplot_i/src
GP_OBJ = $(GP)/gnuplot_i.o

# --- Primary test target ---

TEST_SRCS = test_minidsp.c test_error.c test_core.c test_effects.c \
            test_generators.c test_spectrum.c test_fir.c test_gcc.c \
            test_biquad.c test_windows.c test_dtmf.c test_spectext.c \
            test_steg.c test_fileio.c test_resample.c test_vad.c

TEST_OBJS = $(TEST_SRCS:.c=.o)

HEADERS = ../include/minidsp.h ../include/biquad.h ../include/fileio.h test_helpers.h

test_minidsp: $(TEST_OBJS) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lfftw3 -lsndfile -lm

$(TEST_OBJS): %.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Convenience target: build and run
.PHONY: test
test: test_minidsp
	./test_minidsp

# --- Legacy demo programs ---

gcc_phat_test: gcc_phat_test.o $(GP_OBJ) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ gcc_phat_test.o $(GP_OBJ) $(LDFLAGS) -L.. -lminidsp -lfftw3 -lm

gcc_phat_test.o: gcc_phat_test.c
	$(CC) $(CFLAGS) -I $(GP) -c $< -o $@

gcc_phat_multi_test: gcc_phat_multi_test.o ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lfftw3 -lm

gcc_phat_multi_test.o: gcc_phat_multi_test.c
	$(CC) $(CFLAGS) -c $< -o $@

testliveio: testliveio.o ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) -lportaudio -lm

testliveio.o: testliveio.c
	$(CC) $(CFLAGS) -c $< -o $@

testliverec: testliverec.o $(GP_OBJ) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ testliverec.o $(GP_OBJ) $(LDFLAGS) -L.. -lminidsp -lportaudio -lm

testliverec.o: testliverec.c
	$(CC) $(CFLAGS) -I $(GP) -c $< -o $@

# --- Housekeeping ---

.PHONY: clean
clean:
	-rm -f *.o test_minidsp gcc_phat_test gcc_phat_multi_test testliveio testliverec
	-rm -rf *.dSYM

default: test_minidsp
