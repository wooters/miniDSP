
# Makefile for miniDSP examples
#
# Build all examples:    make
# Run and plot:          make plot
# Clean:                 make clean

include ../config.mk

CFLAGS += -O2 -g -I..

EXAMPLES = magnitude_spectrum power_spectral_density spectrogram phase_spectrum sine_wave white_noise impulse chirp_wave

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../libminidsp.a
../libminidsp.a:
	$(MAKE) -C .. libminidsp.a

default: $(EXAMPLES)

$(EXAMPLES): %: %.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

# Run the examples and generate the interactive HTML plots
.PHONY: plot
plot: $(EXAMPLES)
	./magnitude_spectrum
	./power_spectral_density
	./spectrogram
	./phase_spectrum
	./sine_wave
	./white_noise
	./impulse
	./chirp_wave

.PHONY: clean
clean:
	-rm -rf $(EXAMPLES) $(addsuffix .dSYM,$(EXAMPLES)) $(addsuffix .csv,$(EXAMPLES)) $(addsuffix .html,$(EXAMPLES))
