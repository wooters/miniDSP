# Simple Effects {#simple-effects}

Three classic audio effects built from short delay lines:

- [Delay/Echo](https://en.wikipedia.org/wiki/Delay_%28audio_effect%29) ([circular buffer](https://en.wikipedia.org/wiki/Circular_buffer) + feedback)
- [Tremolo](https://en.wikipedia.org/wiki/Tremolo) (low-frequency [amplitude modulation](https://en.wikipedia.org/wiki/Amplitude_modulation))
- [Comb-filter](https://en.wikipedia.org/wiki/Comb_filter) [reverb](https://en.wikipedia.org/wiki/Reverberation) (feedback comb)

These are simple enough to study sample-by-sample, but they are also the
building blocks of many larger audio effects.

---

## Delay line / echo

Delay/echo reads an older sample from a circular delay line and mixes it
with the current input:

\f[
s[n] = x[n] + feedback \cdot s[n-D]
\f]

The output mixes dry input with the delayed state:

\f[
y[n] = dry \cdot x[n] + wet \cdot s[n-D]
\f]

where \f$D\f$ is the delay in samples and \f$|feedback| < 1\f$.

**Reading the formula in C:**

```c
// x[n] -> in[n], s[n-D] -> d (delay[idx]), y[n] -> out[n], D -> delay_samples
double d = delay[idx];
out[n] = dry * in[n] + wet * d;
delay[idx] = in[n] + feedback * d;
idx = (idx + 1) % delay_samples;
```

**API:**

```c
void MD_delay_echo(const double *in, double *out, unsigned N,
                   unsigned delay_samples, double feedback,
                   double dry, double wet);
```

**Quick example:**

\snippet simple_effects.c delay-echo

**Audio (before/after):**

\htmlonly
<p><strong>Before (dry click train)</strong></p>
<audio controls style="margin: 0.5em 0;">
  <source src="effect_delay_before.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
<p><strong>After (echo: delay=11025, feedback=0.45, dry=1.0, wet=0.6)</strong></p>
<audio controls style="margin: 0.5em 0 1em 0;">
  <source src="effect_delay_after.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Spectrograms (before/after):**

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="effect_delay_before_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="effect_delay_after_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

---

## Tremolo

Tremolo is amplitude modulation by a low-frequency oscillator (LFO).
The gain is:

\f[
g[n] = (1-depth) + depth \cdot \frac{1 + \sin(2\pi f_{LFO}n/f_s)}{2}
\f]

and output is:

\f[
y[n] = g[n] \cdot x[n]
\f]

So gain moves between \f$1-depth\f$ and \f$1\f$.

**Reading the formula in C:**

```c
// f_LFO -> rate_hz, fs -> sample_rate, g[n] -> gain, x[n] -> in[n], y[n] -> out[n]
double lfo = 0.5 * (1.0 + sin(2.0 * M_PI * rate_hz * n / sample_rate));
double gain = (1.0 - depth) + depth * lfo;
out[n] = in[n] * gain;
```

**API:**

```c
void MD_tremolo(const double *in, double *out, unsigned N,
                double rate_hz, double depth, double sample_rate);
```

**Quick example:**

\snippet simple_effects.c tremolo

**Audio (before/after):**

\htmlonly
<p><strong>Before (dry sine, 220 Hz)</strong></p>
<audio controls style="margin: 0.5em 0;">
  <source src="effect_tremolo_before.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
<p><strong>After (tremolo: rate=5.0 Hz, depth=0.8)</strong></p>
<audio controls style="margin: 0.5em 0 1em 0;">
  <source src="effect_tremolo_after.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Spectrograms (before/after):**

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="effect_tremolo_before_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="effect_tremolo_after_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

---

## Comb-filter reverb

A feedback comb filter reuses a delayed copy of its own output:

\f[
c[n] = x[n] + feedback \cdot c[n-D]
\f]

Then we blend dry input and comb output:

\f[
y[n] = dry \cdot x[n] + wet \cdot c[n]
\f]

The repeated, closely spaced echoes create a reverb-like resonant tail.

**Reading the formula in C:**

```c
// c[n-D] -> delayed (comb[idx]), c[n] -> c, x[n] -> in[n], y[n] -> out[n], D -> delay_samples
double delayed = comb[idx];
double c = in[n] + feedback * delayed;
comb[idx] = c;
out[n] = dry * in[n] + wet * c;
idx = (idx + 1) % delay_samples;
```

**API:**

```c
void MD_comb_reverb(const double *in, double *out, unsigned N,
                    unsigned delay_samples, double feedback,
                    double dry, double wet);
```

**Quick example:**

\snippet simple_effects.c comb-reverb

**Audio (before/after):**

\htmlonly
<p><strong>Before (dry decaying tone burst)</strong></p>
<audio controls style="margin: 0.5em 0;">
  <source src="effect_comb_before.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
<p><strong>After (comb reverb: delay=1323, feedback=0.75, dry=0.7, wet=0.6)</strong></p>
<audio controls style="margin: 0.5em 0 1em 0;">
  <source src="effect_comb_after.wav" type="audio/wav">
  <em>Your browser does not support the audio element.</em>
</audio>
\endhtmlonly

**Spectrograms (before/after):**

\htmlonly
<div style="display:flex;gap:0.75rem;margin:1em 0;flex-wrap:wrap;">
  <iframe src="effect_comb_before_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
  <iframe src="effect_comb_after_spectrogram.html" style="flex:1;min-width:280px;height:380px;border:1px solid #ddd;border-radius:4px;" frameborder="0"></iframe>
</div>
\endhtmlonly

---

## Verification tips

- Delay echo impulse input should repeat every `delay_samples`, decaying by `feedback^m`.
- Tremolo with `depth = 0` should exactly match the input.
- Comb reverb with `dry = 0` and impulse input should produce a geometric series at delay multiples.

## API reference

- MD_delay_echo() -- circular-buffer echo with feedback
- MD_tremolo() -- sinusoidal amplitude modulation
- MD_comb_reverb() -- feedback comb filter with dry/wet mix
