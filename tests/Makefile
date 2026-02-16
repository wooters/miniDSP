
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

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -pedantic -g -I..

# Paths to optional dependencies (legacy demos)
GP = $(HOME)/src/gnuplot_i/src
GP_OBJ = $(GP)/gnuplot_i.o

# --- Primary test target ---

test_minidsp: test_minidsp.o ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ -lfftw3 -lm

test_minidsp.o: test_minidsp.c ../minidsp.h ../biquad.h
	$(CC) $(CFLAGS) -c $< -o $@

# Convenience target: build and run
.PHONY: test
test: test_minidsp
	./test_minidsp

# --- Legacy demo programs ---

gcc_phat_test: gcc_phat_test.o $(GP_OBJ) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ gcc_phat_test.o $(GP_OBJ) -L.. -lminidsp -lfftw3 -lm

gcc_phat_test.o: gcc_phat_test.c
	$(CC) $(CFLAGS) -I $(GP) -c $< -o $@

gcc_phat_multi_test: gcc_phat_multi_test.o ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ -lfftw3 -lm

gcc_phat_multi_test.o: gcc_phat_multi_test.c
	$(CC) $(CFLAGS) -c $< -o $@

testliveio: testliveio.o ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $^ -lportaudio -lm

testliveio.o: testliveio.c
	$(CC) $(CFLAGS) -c $< -o $@

testliverec: testliverec.o $(GP_OBJ) ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ testliverec.o $(GP_OBJ) -L.. -lminidsp -lportaudio -lm

testliverec.o: testliverec.c
	$(CC) $(CFLAGS) -I $(GP) -c $< -o $@

# --- Housekeeping ---

.PHONY: clean
clean:
	-rm -f *.o test_minidsp gcc_phat_test gcc_phat_multi_test testliveio testliverec

default: test_minidsp
