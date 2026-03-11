
# Makefile for miniDSP examples
#
# Build all examples:    make
# Run and plot:          make plot
# Clean:                 make clean

include ../config.mk

CFLAGS += -O2 -g -I../include

EXAMPLES = magnitude_spectrum power_spectral_density spectrogram phase_spectrum sine_wave white_noise impulse chirp_wave square_sawtooth basic_operations window_functions fir_convolution simple_effects pitch_detection mel_mfcc

# Examples that need libsndfile (for WAV file I/O)
SNDFILE_EXAMPLES = dtmf_detector spectrogram_text shepard_tone

ALL_EXAMPLES = $(EXAMPLES) $(SNDFILE_EXAMPLES)

# Always delegate to parent Makefile to check if library needs rebuilding
.PHONY: ../libminidsp.a
../libminidsp.a:
	$(MAKE) -C .. libminidsp.a

default: $(ALL_EXAMPLES)

$(EXAMPLES): %: %.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

$(SNDFILE_EXAMPLES): %: %.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

# Docs-only utilities (not user-facing examples — built by `make docs` from root)
gen_audio_samples: gen_audio_samples.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lsndfile -lm

gen_signal_plots: gen_signal_plots.c ../libminidsp.a
	$(CC) $(CFLAGS) -o $@ $< -L.. -lminidsp $(LDFLAGS) -lfftw3 -lm

# Run the examples and generate the interactive HTML plots
.PHONY: plot
plot: $(ALL_EXAMPLES)
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
	./pitch_detection
	./mel_mfcc
	./dtmf_detector
	./spectrogram_text
	./shepard_tone

.PHONY: clean
clean:
	-rm -rf $(ALL_EXAMPLES) $(addsuffix .dSYM,$(ALL_EXAMPLES)) $(addsuffix .csv,$(ALL_EXAMPLES)) $(addsuffix .html,$(ALL_EXAMPLES))
	-rm -f spectrogram_text.png spectrogram_text.wav shepard_tone.wav
	-rm -rf gen_audio_samples gen_audio_samples.dSYM
	-rm -rf gen_signal_plots gen_signal_plots.dSYM
