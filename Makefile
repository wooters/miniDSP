
# miniDSP Library Makefile
#
# Build the static library:    make
# Build and run tests:         make test
# Generate documentation:      make docs
# Clean everything:            make clean

include config.mk

CFLAGS += -O2 -g -Iinclude

MD_SRCS = src/minidsp.c src/biquad.c src/liveio.c src/fileio.c
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

# Install git hooks (pre-push runs tests before pushing to main)
.PHONY: install-hooks
install-hooks:
	ln -sf ../../scripts/pre-push .git/hooks/pre-push
	@echo "Git hooks installed."

# Generate HTML documentation with Doxygen
.PHONY: docs
docs: libminidsp.a
	$(MAKE) -C examples gen_audio_samples
	$(MAKE) -C examples gen_signal_plots
	mkdir -p guides/audio guides/plots
	examples/gen_audio_samples
	examples/gen_signal_plots
	doxygen Doxyfile

.PHONY: clean
clean:
	-rm -f src/*.o libminidsp.a
	-rm -rf docs guides/audio guides/plots
	$(MAKE) -C tests clean
	$(MAKE) -C examples clean
