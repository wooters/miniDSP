# TODO -- DSP Feature Gap Checklist

Features that would round out miniDSP for students learning audio DSP,
ranked roughly by educational impact.

## FFT / Spectrum Analysis

- [x] **Magnitude spectrum** -- compute |X(k)| from a real signal; the first thing students plot after learning the DFT
- [ ] **Power spectral density (PSD)** -- power per frequency bin (|X(k)|^2 / N); essential for noise analysis
- [ ] **Spectrogram / STFT** -- sliding-window FFT producing a time-frequency matrix; the workhorse of audio visualization
- [ ] **Phase spectrum** -- arg(X(k)); needed for phase-vocoder effects and understanding group delay

## Signal Generators

- [ ] **Sine generator** -- constant-frequency tone at a given amplitude and sample rate; the "hello world" of DSP
- [ ] **White noise generator** -- uniformly or normally distributed random samples; used to test filters and measure impulse responses
- [ ] **Impulse generator** -- single unit-sample spike; reveals a system's impulse response directly
- [ ] **Chirp (swept sine)** -- frequency sweeps linearly or logarithmically over time; great for testing filter magnitude response
- [ ] **Square / sawtooth waves** -- classic non-sinusoidal waveforms; demonstrate harmonics and Gibbs phenomenon

## FIR Filters / Convolution

- [ ] **Time-domain convolution** -- direct sum-of-products implementation; teaches the convolution theorem before FFT shortcuts
- [ ] **Moving-average filter** -- simplest FIR low-pass; good first filter for students to implement and analyse
- [ ] **General FIR filter** -- apply arbitrary coefficient arrays; pairs with window-design method
- [ ] **FFT-based fast convolution** -- overlap-add or overlap-save; shows why FFT matters for long filters

## Window Functions

Hanning is already implemented (`MD_hanning_window`). Add the other common ones
so students can compare sidelobe trade-offs.

- [ ] **Hamming window** -- slightly lower first sidelobe than Hanning; popular in speech processing
- [ ] **Blackman window** -- much lower sidelobes at the cost of a wider main lobe
- [ ] **Rectangular window** -- trivial (all ones) but useful as a baseline reference

## Basic Signal Operations

- [ ] **RMS (root mean square)** -- standard loudness measure; more perceptually meaningful than peak amplitude
- [ ] **Zero-crossing rate** -- simple proxy for pitch/noisiness; widely used in speech/music classification
- [ ] **Autocorrelation** -- similarity of a signal with a delayed copy of itself; foundation of pitch detection
- [ ] **Peak detection** -- find local maxima above a threshold; used in onset detection and pitch tracking
- [ ] **Signal mixing** -- weighted sum of two signals; needed for any multi-source demo

## Simple Effects

- [ ] **Delay line / echo** -- circular-buffer delay with feedback; the building block of most audio effects
- [ ] **Tremolo** -- amplitude modulation by a low-frequency oscillator; simplest audible effect to implement
- [ ] **Comb-filter reverb** -- feedback comb filter; introduces students to IIR reverb design (Schroeder reverberator)

## Pitch Detection

- [ ] **F0 estimation (autocorrelation)** -- find the fundamental frequency by locating the autocorrelation peak; classic speech/music analysis
- [ ] **F0 estimation (FFT-based)** -- peak-pick the magnitude spectrum; simpler but less robust than autocorrelation

## Mel-scale / MFCCs

- [ ] **Mel filterbank** -- triangular filters spaced on the mel scale; bridges DSP and machine-learning feature extraction
- [ ] **MFCCs (Mel-frequency cepstral coefficients)** -- DCT of log mel energies; the standard front-end for speech and audio ML
