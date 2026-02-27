# Shepard Tone Generator {#shepard-tone}

A [Shepard tone](https://en.wikipedia.org/wiki/Shepard_tone) is the
auditory equivalent of an [M.C. Escher](https://en.wikipedia.org/wiki/M._C._Escher)
staircase: a sound that seems to rise (or fall) in pitch forever without
ever actually leaving its frequency range.  The effect was first described
by cognitive scientist [Roger Shepard](https://en.wikipedia.org/wiki/Roger_Shepard)
in 1964.

miniDSP provides a single-call generator in `src/minidsp_generators.c`,
demonstrated in `examples/shepard_tone.c`.

Build and run from the repository root:

```sh
make -C examples shepard_tone
cd examples && ./shepard_tone
```

---

## How it works

The illusion rests on two ideas:

1. **Octave equivalence** — the human ear perceives tones one octave apart
   as the "same note" at a different pitch height.
2. **Spectral envelope** — a fixed [Gaussian](https://en.wikipedia.org/wiki/Gaussian_function)
   bell curve in log-frequency space controls how loud each tone is.
   Tones near the centre are loud; tones at the edges are nearly silent.

Several sine waves are sounded simultaneously, each separated by one octave.
All of them glide upward (or downward) in pitch at the same rate.  As a tone
approaches the upper edge of the Gaussian, it fades to silence.  Meanwhile,
a new tone enters at the lower edge, fading in.  Because the only thing the
ear can latch onto (the **loudest** tones) are always in the middle and always
going up, the sound appears to ascend indefinitely.

This diagram shows the principle for a **rising** Shepard tone:

```
Amplitude
  |                ╭───╮
  |              ╱       ╲          ← Gaussian spectral envelope
  |            ╱     ●     ╲          (fixed in log-frequency)
  |          ╱    ●     ●    ╲
  |        ╱   ●           ●   ╲
  |      ╱  ●                 ●  ╲
  |    ╱ ●                       ● ╲
  +--●─────────────────────────────●──→ log₂(frequency)
     ↑                               ↑
  low edge                       high edge
  (fading in)                  (fading out)

  ● = individual octave-spaced tones, all gliding →
```

---

## Signal model

At time \f$t = n / f_s\f$, the output sample is:

\f[
x[n] \;=\; A_\text{norm}\,\sum_k\;
  \underbrace{
    \exp\!\Bigl(-\frac{d_k(t)^2}{2\sigma^2}\Bigr)
  }_{\text{Gaussian envelope}}
  \;\sin\!\bigl(\varphi_k(n)\bigr)
\f]

where:

| Symbol | Meaning |
|--------|---------|
| \f$k\f$ | Layer index (one per octave) |
| \f$d_k(t) = k - c + R\,t\f$ | Octave distance from the Gaussian centre at time \f$t\f$ |
| \f$c = (L-1)/2\f$ | Centre of the layer range |
| \f$\sigma = L/4\f$ | Gaussian width in octaves |
| \f$R\f$ | Glissando rate (`rate_octaves_per_sec`): positive = rising |
| \f$L\f$ | Number of audible octave layers (`num_octaves`) |
| \f$f_k(t) = f_\text{base} \cdot 2^{d_k(t)}\f$ | Instantaneous frequency of layer \f$k\f$ |
| \f$\varphi_k(n)\f$ | Phase accumulated sample-by-sample from \f$f_k\f$ |
| \f$A_\text{norm}\f$ | Normalisation factor so peak amplitude equals the requested `amplitude` |

The Gaussian is fixed in log-frequency space while the tones glide
through it.  Only layers within \f$\pm 5\sigma\f$ are computed;
layers whose frequency exceeds the [Nyquist frequency](https://en.wikipedia.org/wiki/Nyquist_frequency)
or falls below 20 Hz are silently skipped.

---

## Reading the formula in C

The core synthesis loop — what `MD_shepard_tone()` does internally:

```c
// R -> rate,  L -> num_octaves,  f_base -> base_freq,  fs -> sample_rate
// c = (L - 1) / 2,  sigma = L / 4
// d_k(t) -> d,  f_k(t) -> freq,  phi_k(n) -> phases[k]
// x[n] -> output[i]

double c     = (double)(num_octaves - 1) / 2.0;
double sigma = (double)num_octaves / 4.0;

for (unsigned i = 0; i < N; i++) {
    double t = (double)i / sample_rate;
    double sample = 0.0;
    for (int k = k_min; k <= k_max; k++) {
        double d    = (double)k - c + rate * t;       // octave distance from centre
        double freq = base_freq * pow(2.0, d);        // instantaneous frequency
        double gauss = exp(-0.5 * d * d / (sigma * sigma));  // Gaussian weight
        phases[k] += 2.0 * M_PI * freq / sample_rate; // accumulate phase
        sample += gauss * sin(phases[k]);              // add weighted sine
    }
    output[i] = sample;  // (later normalised to peak amplitude)
}
```

---

## Parameters and their effect

### Glissando rate (rate_octaves_per_sec)

Controls how fast the tones rise or fall.

| Rate | Effect |
|------|--------|
| 0.0 | Static chord — no movement, just octave-spaced sines |
| 0.25 | Slow, dreamy ascent (4 seconds per octave) |
| 0.5 | Moderate rise (default) — 2 seconds per octave |
| 1.0 | Fast rise — 1 second per octave |
| −0.5 | Moderate descent — the "falling" Shepard tone |

**Listen** — rising at 0.5 oct/s (5 seconds):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="shepard_rising.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Listen** — falling at 0.5 oct/s (5 seconds):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="shepard_falling.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Listen** — static chord (5 seconds):

\htmlonly
<audio controls style="margin: 0.5em 0;">
  <source src="shepard_static.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

### Rising spectrogram

The spectrogram below shows the characteristic pattern of a rising Shepard
tone: parallel diagonal lines (one per octave layer) sweeping upward through
the Gaussian bell curve.  Tones fade in at the bottom and fade out at the top.

\htmlonly
<iframe src="shepard_rising_spectrogram.html" style="width:100%;height:420px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

### Falling spectrogram

The falling variant mirrors the rising pattern — diagonal lines sweep
**downward**.

\htmlonly
<iframe src="shepard_falling_spectrogram.html" style="width:100%;height:420px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
\endhtmlonly

---

### Number of octaves (num_octaves)

Controls how many simultaneous octave layers are present and the width
of the Gaussian envelope (\f$\sigma = L/4\f$).

| Value | Typical use |
|-------|-------------|
| 4–6 | Narrow bell — prominent entry/exit of tones; more "organ-like" |
| 8 | Default — smooth, balanced illusion |
| 10–12 | Wide bell — very gradual fading; ethereal, diffuse texture |

More octaves means more layers span the audible range at any instant,
making the transitions smoother at the expense of a busier spectrum.

### Base frequency (base_freq)

The centre of the Gaussian bell curve.  Tones above this frequency are
treated the same as those below it (the Gaussian is symmetric in
log-frequency space).  Typical values: 200–600 Hz.

---

## API

```c
void MD_shepard_tone(double *output, unsigned N, double amplitude,
                     double base_freq, double sample_rate,
                     double rate_octaves_per_sec, unsigned num_octaves);
```

**Parameters:**

| Parameter                | Description |
|--------------------------|-------------|
| `output`                 | Caller-allocated buffer for the synthesised audio. |
| `N`                      | Number of samples to generate.  Must be > 0. |
| `amplitude`              | Peak amplitude of the output signal. |
| `base_freq`              | Centre frequency of the Gaussian envelope in Hz (e.g. 440). |
| `sample_rate`            | Sample rate in Hz.  Must be > 0. |
| `rate_octaves_per_sec`   | Glissando rate: positive = rising, negative = falling, 0 = static. |
| `num_octaves`            | Number of audible octave layers (Gaussian width).  Must be > 0. |

---

## Quick example

```c
#include "minidsp.h"
#include <stdlib.h>

// 5 seconds of endlessly rising Shepard tone at 44.1 kHz
unsigned N = 5 * 44100;
double *sig = malloc(N * sizeof(double));
MD_shepard_tone(sig, N, 0.8, 440.0, 44100.0, 0.5, 8);
// sig[] now sounds like it rises forever
free(sig);
```

---

## Example program

The example `examples/shepard_tone.c` generates a WAV file and an interactive
HTML spectrogram.

**Usage:**

```sh
./shepard_tone [--rising | --falling | --static]
               [--rate OCTAVES_PER_SEC]
               [--octaves NUM]
               [--base FREQ_HZ]
               [--duration SEC]
```

Default: rising at 0.5 oct/s, 440 Hz base, 8 octaves, 5 seconds.

**Generate and listen:**

\snippet shepard_tone.c generate-signal

**Build and run:**

```sh
make -C examples shepard_tone
cd examples && ./shepard_tone
open shepard_tone.html     # interactive spectrogram
```

---

## Why it works — the psychoacoustics

The Shepard tone exploits a fundamental ambiguity in pitch perception.
Pitch has two dimensions:

- **[Pitch chroma](https://en.wikipedia.org/wiki/Pitch_class)** — which note
  it is (C, D, E, …), determined by the position within the octave.
- **Pitch height** — how high or low it sounds overall.

The Gaussian envelope removes pitch-height cues: there is no single
"highest" or "lowest" tone to anchor the percept.  All the listener
hears is the chroma — and the chroma is always going up.

The effect is even more striking when a rising Shepard tone is followed
by a falling one.  Despite the falling version being physically the
mirror image, many listeners perceive *both* as rising — a dramatic
demonstration of how expectation shapes perception.

[Jean-Claude Risset](https://en.wikipedia.org/wiki/Jean-Claude_Risset) later
extended the idea to **continuous glissando** (the Shepard–Risset glissando),
which is exactly what `MD_shepard_tone()` implements: instead of discrete
steps, the tones slide smoothly.

---

## Further reading

- [Roger N. Shepard, "Circularity in Judgments of Relative Pitch" (1964)](https://doi.org/10.1121/1.1908239) — the original paper.
- [Shepard tone — Wikipedia](https://en.wikipedia.org/wiki/Shepard_tone)
- [Shepard–Risset glissando — Wikipedia](https://en.wikipedia.org/wiki/Shepard_tone#Shepard%E2%80%93Risset_glissando)
- Diana Deutsch, [Musical Illusions and Phantom Words](https://deutsch.ucsd.edu/psychology/pages.php?i=201) — demonstrations and further reading on auditory illusions.
