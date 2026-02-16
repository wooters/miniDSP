
# miniDSP Library Makefile
#
# Build the static library:    make
# Build and run tests:         make test
# Clean everything:            make clean

CC = gcc
CFLAGS = -std=c11 -Wall -Wextra -pedantic -O2 -g

MD_SRCS = minidsp.c biquad.c liveio.c fileio.c
MD_OBJS = $(MD_SRCS:.c=.o)

default: libminidsp.a

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

libminidsp.a: $(MD_OBJS)
	ar rcs $@ $(MD_OBJS)

# Build and run the test suite
.PHONY: test
test: libminidsp.a
	$(MAKE) -C tests test

.PHONY: clean
clean:
	-rm -f *.o libminidsp.a
	$(MAKE) -C tests clean
