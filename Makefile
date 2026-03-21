
# miniDSP Library Makefile
#
# Build the static library:    make
# Build and run tests:         make test
# Generate documentation:      make docs
# Clean everything:            make clean

include config.mk

VERSION := $(shell cat VERSION)
VERSION_MAJOR := $(word 1,$(subst ., ,$(VERSION)))
VERSION_MINOR := $(word 2,$(subst ., ,$(VERSION)))
VERSION_PATCH := $(word 3,$(subst ., ,$(VERSION)))

CFLAGS += -O2 -g -Iinclude
CFLAGS += -DMINIDSP_VERSION_MAJOR=$(VERSION_MAJOR) \
          -DMINIDSP_VERSION_MINOR=$(VERSION_MINOR) \
          -DMINIDSP_VERSION_PATCH=$(VERSION_PATCH) \
          -DMINIDSP_VERSION='"$(VERSION)"'

MD_SRCS = src/minidsp_error.c src/minidsp_core.c src/minidsp_generators.c \
          src/minidsp_spectrum.c src/minidsp_fir.c src/minidsp_dtmf.c \
          src/minidsp_spectext.c src/minidsp_steg.c src/minidsp_gcc.c \
          src/minidsp_resample.c src/minidsp_vad.c src/biquad.c \
          src/liveio.c src/fileio.c
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

# Build tools
.PHONY: tools
tools: libminidsp.a
	$(MAKE) -C tools/mel_viz
	$(MAKE) -C tools/audio_steg
	$(MAKE) -C tools/resample

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
	$(MAKE) -C examples vad
	mkdir -p guides/audio guides/plots
	examples/gen_audio_samples
	examples/gen_signal_plots
	cd examples && ./vad
	doxygen Doxyfile
	python3 scripts/gen_llms_txt.py

PREFIX ?= /usr/local
HEADERS = $(wildcard include/*.h)

.PHONY: install
install: libminidsp.a
	install -d $(PREFIX)/lib $(PREFIX)/include
	install -m 644 libminidsp.a $(PREFIX)/lib/
	install -m 644 $(HEADERS) $(PREFIX)/include/

.PHONY: uninstall
uninstall:
	rm -f $(PREFIX)/lib/libminidsp.a
	$(foreach h,$(notdir $(HEADERS)),rm -f $(PREFIX)/include/$(h);)

.PHONY: clean
clean:
	-rm -f src/*.o libminidsp.a
	-rm -rf docs/html docs/xml guides/audio guides/plots
	$(MAKE) -C tests clean
	$(MAKE) -C examples clean
	$(MAKE) -C tools/mel_viz clean
	$(MAKE) -C tools/audio_steg clean
	$(MAKE) -C tools/resample clean
