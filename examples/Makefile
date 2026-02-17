
# Makefile for miniDSP examples
#
# Build all examples:    make
# Run and plot:          make plot
# Clean:                 make clean

CC = gcc
CFLAGS = -std=c23 -Wall -Wextra -pedantic -O2 -g -I..

# Homebrew on Apple Silicon installs to /opt/homebrew; on Intel to /usr/local.
BREW_PREFIX := $(shell brew --prefix 2>/dev/null)
ifneq ($(BREW_PREFIX),)
  CFLAGS += -I$(BREW_PREFIX)/include
  LDFLAGS += -L$(BREW_PREFIX)/lib
endif

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../libminidsp.a
../libminidsp.a:
	$(MAKE) -C .. libminidsp.a

default: magnitude_spectrum power_spectral_density

magnitude_spectrum: magnitude_spectrum.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

power_spectral_density: power_spectral_density.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

# Run the examples and generate the interactive HTML plots
.PHONY: plot
plot: magnitude_spectrum power_spectral_density
	./magnitude_spectrum
	./power_spectral_density

.PHONY: clean
clean:
	-rm -rf magnitude_spectrum magnitude_spectrum.dSYM magnitude_spectrum.csv magnitude_spectrum.html
	-rm -rf power_spectral_density power_spectral_density.dSYM power_spectral_density.csv power_spectral_density.html
