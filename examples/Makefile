
# Makefile for miniDSP examples
#
# Build all examples:    make
# Run and plot:          make plot
# Clean:                 make clean

include ../config.mk

CFLAGS += -O2 -g -I../include

EXAMPLES = magnitude_spectrum power_spectral_density spectrogram phase_spectrum sine_wave white_noise impulse chirp_wave square_sawtooth basic_operations window_functions fir_convolution simple_effects

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../libminidsp.a
../libminidsp.a:
	$(MAKE) -C .. libminidsp.a

default: $(EXAMPLES)

$(EXAMPLES): %: %.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

# Docs-only utilities (not user-facing examples — built by `make docs` from root)
gen_audio_samples: gen_audio_samples.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

gen_signal_plots: gen_signal_plots.c ../libminidsp.a
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
	./square_sawtooth
	./basic_operations
	./window_functions
	./fir_convolution
	./simple_effects

.PHONY: clean
clean:
	-rm -rf $(EXAMPLES) $(addsuffix .dSYM,$(EXAMPLES)) $(addsuffix .csv,$(EXAMPLES)) $(addsuffix .html,$(EXAMPLES))
	-rm -rf gen_audio_samples gen_audio_samples.dSYM
	-rm -rf gen_signal_plots gen_signal_plots.dSYM
