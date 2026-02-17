
# miniDSP Library Makefile
#
# Build the static library:    make
# Build and run tests:         make test
# Generate documentation:      make docs
# Clean everything:            make clean

CC = gcc
CFLAGS = -std=c23 -Wall -Wextra -pedantic -O2 -g

# Homebrew on Apple Silicon installs to /opt/homebrew; on Intel to /usr/local.
# If Homebrew isn't present, these lines are skipped and system paths are used.
BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
  CFLAGS += -I$(BREW_PREFIX)/include
endif

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

# Build the examples
.PHONY: examples
examples: libminidsp.a
	$(MAKE) -C examples

# Build and test inside an Ubuntu 24.04 container (requires macOS container CLI)
CONTAINER_TAG = minidsp-test

.PHONY: container-test
container-test:
	container build --tag $(CONTAINER_TAG) --file Dockerfile .
	container run --rm $(CONTAINER_TAG)

# Generate HTML documentation with Doxygen
.PHONY: docs
docs:
	doxygen Doxyfile

.PHONY: clean
clean:
	-rm -f *.o libminidsp.a
	-rm -rf docs
	$(MAKE) -C tests clean
	$(MAKE) -C examples clean
