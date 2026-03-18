/**
 * @file mel_viz.c
 * @brief Mel-spectrum audio visualizer — generates browser-based radial animation.
 *
 * Reads a WAV file, computes per-frame mel energies using MD_mel_energies(),
 * and assembles an output folder with mel data + web assets for browser playback.
 *
 * Usage:
 *   mel_viz <input.wav> [options]
 *     -o <dir>         Output directory (default: mel_viz_out/)
 *     --mels <n>       Number of mel bands (default: 24)
 *     --fft-size <n>   FFT window size (default: 2048)
 *     --fps <n>        Frames per second (default: 30)
 *     --min-freq <hz>  Low frequency bound (default: 40)
 *     --max-freq <hz>  High frequency bound (default: 16000)
 *     --width <px>     Canvas width (default: 1080)
 *     --height <px>    Canvas height (default: 1080)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>
#include <sndfile.h>
#include "minidsp.h"

/* -------------------------------------------------------------------
 * Configuration (populated from CLI arguments)
 * ---------------------------------------------------------------- */

typedef struct {
    const char *input_path;
    const char *output_dir;
    unsigned    num_mels;
    unsigned    fft_size;
    unsigned    fps;
    double      min_freq;
    double      max_freq;
    unsigned    width;
    unsigned    height;
} config_t;

static void config_defaults(config_t *cfg)
{
    cfg->input_path = NULL;
    cfg->output_dir = "mel_viz_out";
    cfg->num_mels   = 24;
    cfg->fft_size   = 2048;
    cfg->fps        = 30;
    cfg->min_freq   = 40.0;
    cfg->max_freq   = 16000.0;
    cfg->width      = 1080;
    cfg->height     = 1080;
}

static void usage(const char *progname)
{
    fprintf(stderr,
        "Usage: %s <input.wav> [options]\n"
        "  -o <dir>         Output directory (default: mel_viz_out/)\n"
        "  --mels <n>       Number of mel bands (default: 24)\n"
        "  --fft-size <n>   FFT window size (default: 2048)\n"
        "  --fps <n>        Frames per second (default: 30)\n"
        "  --min-freq <hz>  Low frequency bound (default: 40)\n"
        "  --max-freq <hz>  High frequency bound (default: 16000)\n"
        "  --width <px>     Canvas width (default: 1080)\n"
        "  --height <px>    Canvas height (default: 1080)\n",
        progname);
}

static int parse_args(int argc, char **argv, config_t *cfg)
{
    config_defaults(cfg);

    for (int i = 1; i < argc; i++) {
        if (argv[i][0] != '-') {
            if (cfg->input_path != NULL) {
                fprintf(stderr, "Error: multiple input files specified\n");
                return -1;
            }
            cfg->input_path = argv[i];
        } else if (strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            cfg->output_dir = argv[++i];
        } else if (strcmp(argv[i], "--mels") == 0 && i + 1 < argc) {
            cfg->num_mels = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--fft-size") == 0 && i + 1 < argc) {
            cfg->fft_size = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--fps") == 0 && i + 1 < argc) {
            cfg->fps = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--min-freq") == 0 && i + 1 < argc) {
            cfg->min_freq = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-freq") == 0 && i + 1 < argc) {
            cfg->max_freq = atof(argv[++i]);
        } else if (strcmp(argv[i], "--width") == 0 && i + 1 < argc) {
            cfg->width = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "--height") == 0 && i + 1 < argc) {
            cfg->height = (unsigned)atoi(argv[++i]);
        } else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            usage(argv[0]);
            exit(0);
        } else {
            fprintf(stderr, "Error: unknown option '%s'\n", argv[i]);
            usage(argv[0]);
            return -1;
        }
    }

    if (cfg->input_path == NULL) {
        fprintf(stderr, "Error: no input WAV file specified\n");
        usage(argv[0]);
        return -1;
    }

    return 0;
}

/* -------------------------------------------------------------------
 * WAV I/O
 * ---------------------------------------------------------------- */

/**
 * Read a WAV file into a mono double buffer.
 * If stereo, averages channels. Returns sample count (0 on error).
 */
static unsigned read_wav(const char *path, double **out_samples,
                         double *out_sample_rate)
{
    SF_INFO info;
    memset(&info, 0, sizeof(info));

    SNDFILE *sf = sf_open(path, SFM_READ, &info);
    if (!sf) {
        fprintf(stderr, "Error: cannot open '%s': %s\n", path, sf_strerror(NULL));
        return 0;
    }

    sf_count_t total = info.frames * info.channels;
    double *raw = malloc((size_t)total * sizeof(double));
    if (!raw) {
        fprintf(stderr, "Error: allocation failed for %lld samples\n",
                (long long)total);
        sf_close(sf);
        return 0;
    }

    sf_count_t read = sf_readf_double(sf, raw, info.frames);
    sf_close(sf);

    if (read != info.frames) {
        fprintf(stderr, "Warning: expected %lld frames, read %lld\n",
                (long long)info.frames, (long long)read);
    }

    /* Downmix to mono if needed */
    double *mono = malloc((size_t)info.frames * sizeof(double));
    if (!mono) {
        fprintf(stderr, "Error: allocation failed\n");
        free(raw);
        return 0;
    }

    if (info.channels == 1) {
        memcpy(mono, raw, (size_t)info.frames * sizeof(double));
    } else {
        for (sf_count_t i = 0; i < info.frames; i++) {
            double sum = 0.0;
            for (int ch = 0; ch < info.channels; ch++) {
                sum += raw[i * info.channels + ch];
            }
            mono[i] = sum / info.channels;
        }
    }

    free(raw);
    *out_samples = mono;
    *out_sample_rate = (double)info.samplerate;
    return (unsigned)info.frames;
}

/* -------------------------------------------------------------------
 * Mel analysis
 * ---------------------------------------------------------------- */

/**
 * Compute per-frame mel energies and bass envelope.
 *
 * Returns heap-allocated mel_frames (num_frames * num_mels) and
 * bass_env (num_frames). Caller must free both.
 */
static unsigned compute_mel_frames(const double *samples, unsigned num_samples,
                                   double sample_rate, const config_t *cfg,
                                   double **out_mel, double **out_bass)
{
    unsigned hop = (unsigned)(sample_rate / cfg->fps);
    if (hop == 0) hop = 1;

    /* Count frames that fit */
    unsigned num_frames = 0;
    if (num_samples >= cfg->fft_size) {
        num_frames = (num_samples - cfg->fft_size) / hop + 1;
    }

    if (num_frames == 0) {
        fprintf(stderr, "Error: audio too short for FFT size %u\n", cfg->fft_size);
        *out_mel = NULL;
        *out_bass = NULL;
        return 0;
    }

    double *mel_frames = malloc((size_t)num_frames * cfg->num_mels * sizeof(double));
    double *bass_env   = malloc((size_t)num_frames * sizeof(double));
    if (!mel_frames || !bass_env) {
        fprintf(stderr, "Error: allocation failed\n");
        free(mel_frames);
        free(bass_env);
        *out_mel = NULL;
        *out_bass = NULL;
        return 0;
    }

    /* Compute mel energies per frame */
    for (unsigned f = 0; f < num_frames; f++) {
        unsigned offset = f * hop;
        MD_mel_energies(samples + offset, cfg->fft_size,
                        sample_rate, cfg->num_mels,
                        cfg->min_freq, cfg->max_freq,
                        mel_frames + (size_t)f * cfg->num_mels);
    }

    /* Compute bass envelope: sum of lowest 6 mel bands + envelope follower */
    unsigned bass_bands = cfg->num_mels < 6 ? cfg->num_mels : 6;
    double decay = 0.85;  /* envelope follower decay per frame */
    double env = 0.0;

    for (unsigned f = 0; f < num_frames; f++) {
        double bass = 0.0;
        for (unsigned b = 0; b < bass_bands; b++) {
            bass += mel_frames[(size_t)f * cfg->num_mels + b];
        }
        env = bass > env ? bass : env * decay;
        bass_env[f] = env;
    }

    *out_mel = mel_frames;
    *out_bass = bass_env;
    return num_frames;
}

/* -------------------------------------------------------------------
 * Linear analysis
 * ---------------------------------------------------------------- */

/**
 * Compute per-frame linear-spaced band energies.
 *
 * Applies a Hann window, computes one-sided PSD, then bins into num_bands
 * uniformly-spaced frequency bands over [min_freq, max_freq].
 *
 * num_frames must match the value returned by compute_mel_frames() for the
 * same config — caller is responsible for passing the correct count.
 */
static void compute_linear_frames(const double *samples, unsigned num_samples,
                                  double sample_rate, const config_t *cfg,
                                  unsigned num_frames, double *linear_out)
{
    unsigned hop = (unsigned)(sample_rate / cfg->fps);
    if (hop == 0) hop = 1;

    unsigned N = cfg->fft_size;
    unsigned num_bins = N / 2 + 1;
    double bin_hz = sample_rate / N;

    /* Precompute Hann window */
    double *win = malloc(N * sizeof(double));
    double *buf = malloc(N * sizeof(double));
    double *psd = malloc(num_bins * sizeof(double));
    assert(win && buf && psd);
    MD_Gen_Hann_Win(win, N);

    for (unsigned f = 0; f < num_frames; f++) {
        unsigned offset = f * hop;
        assert(offset + N <= num_samples);

        /* Apply Hann window */
        for (unsigned i = 0; i < N; i++)
            buf[i] = samples[offset + i] * win[i];

        /* Compute PSD: |X(k)|^2 / N */
        MD_power_spectral_density(buf, N, psd);

        /* Bin PSD into num_bands uniform Hz bands */
        double *row = linear_out + (size_t)f * cfg->num_mels;
        for (unsigned b = 0; b < cfg->num_mels; b++) {
            double band_lo = cfg->min_freq + (cfg->max_freq - cfg->min_freq) * b / cfg->num_mels;
            double band_hi = cfg->min_freq + (cfg->max_freq - cfg->min_freq) * (b + 1) / cfg->num_mels;
            double sum = 0.0;
            for (unsigned k = 0; k < num_bins; k++) {
                double freq = k * bin_hz;
                if (freq >= band_lo && freq < band_hi)
                    sum += psd[k];
            }
            row[b] = sum;
        }
    }

    free(win);
    free(buf);
    free(psd);
}

/* -------------------------------------------------------------------
 * Output assembly
 * ---------------------------------------------------------------- */

static int mkdir_p(const char *path)
{
    struct stat st;
    if (stat(path, &st) == 0) return 0;
    return mkdir(path, 0755);
}

static int copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    if (!in) {
        fprintf(stderr, "Error: cannot open '%s' for reading\n", src);
        return -1;
    }

    FILE *out = fopen(dst, "wb");
    if (!out) {
        fprintf(stderr, "Error: cannot open '%s' for writing\n", dst);
        fclose(in);
        return -1;
    }

    char buf[8192];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), in)) > 0) {
        if (fwrite(buf, 1, n, out) != n) {
            fprintf(stderr, "Error: write failed to '%s'\n", dst);
            fclose(in);
            fclose(out);
            return -1;
        }
    }

    fclose(in);
    fclose(out);
    return 0;
}

static int write_data_js(const char *path, const config_t *cfg,
                         double sample_rate, unsigned num_frames,
                         const double *mel_frames, const double *linear_frames,
                         const double *bass_env)
{
    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot open '%s' for writing\n", path);
        return -1;
    }

    fprintf(fp, "// Generated by mel_viz\n");
    fprintf(fp, "const MEL_VIZ_DATA = {\n");
    fprintf(fp, "  sampleRate: %.0f,\n", sample_rate);
    fprintf(fp, "  fps: %u,\n", cfg->fps);
    fprintf(fp, "  numFrames: %u,\n", num_frames);
    fprintf(fp, "  numBands: %u,\n", cfg->num_mels);
    fprintf(fp, "  width: %u,\n", cfg->width);
    fprintf(fp, "  height: %u,\n", cfg->height);
    fprintf(fp, "  audioFile: \"audio.wav\",\n");

    /* Write mel frames as flat array */
    fprintf(fp, "  frames: [\n    ");
    size_t total = (size_t)num_frames * cfg->num_mels;
    for (size_t i = 0; i < total; i++) {
        fprintf(fp, "%.6g", mel_frames[i]);
        if (i + 1 < total) {
            fprintf(fp, ",");
            if ((i + 1) % 12 == 0) fprintf(fp, "\n    ");
        }
    }
    fprintf(fp, "\n  ],\n");

    /* Write linear frames as flat array */
    fprintf(fp, "  linearFrames: [\n    ");
    for (size_t i = 0; i < total; i++) {
        fprintf(fp, "%.6g", linear_frames[i]);
        if (i + 1 < total) {
            fprintf(fp, ",");
            if ((i + 1) % 12 == 0) fprintf(fp, "\n    ");
        }
    }
    fprintf(fp, "\n  ],\n");

    /* Write bass envelope */
    fprintf(fp, "  bassEnvelope: [\n    ");
    for (unsigned f = 0; f < num_frames; f++) {
        fprintf(fp, "%.6g", bass_env[f]);
        if (f + 1 < num_frames) {
            fprintf(fp, ",");
            if ((f + 1) % 16 == 0) fprintf(fp, "\n    ");
        }
    }
    fprintf(fp, "\n  ]\n");

    fprintf(fp, "};\n");
    fclose(fp);
    return 0;
}

/**
 * Find the directory containing the mel_viz executable.
 * This is used to locate the web/ assets relative to the binary.
 */
static void get_exe_dir(const char *argv0, char *dir, size_t dir_size)
{
    /* Try to find the last path separator */
    const char *last_sep = strrchr(argv0, '/');
    if (last_sep) {
        size_t len = (size_t)(last_sep - argv0);
        if (len >= dir_size) len = dir_size - 1;
        memcpy(dir, argv0, len);
        dir[len] = '\0';
    } else {
        /* No path separator — binary is in current directory */
        dir[0] = '.';
        dir[1] = '\0';
    }
}

static int assemble_output(const char *argv0, const config_t *cfg,
                           double sample_rate, unsigned num_frames,
                           const double *mel_frames, const double *linear_frames,
                           const double *bass_env)
{
    /* Create output directory */
    if (mkdir_p(cfg->output_dir) != 0) {
        fprintf(stderr, "Error: cannot create directory '%s'\n", cfg->output_dir);
        return -1;
    }

    /* Find web assets relative to executable */
    char exe_dir[4096];
    get_exe_dir(argv0, exe_dir, sizeof(exe_dir));

    char src_path[4096];
    char dst_path[4096];

    /* Copy web assets */
    const char *web_files[] = {
        "index.html", "renderer.js", "audio-provider.js",
        "palettes.js", "controls.js", "exporter.js", "style.css", NULL
    };

    for (int i = 0; web_files[i]; i++) {
        snprintf(src_path, sizeof(src_path), "%s/web/%s", exe_dir, web_files[i]);
        snprintf(dst_path, sizeof(dst_path), "%s/%s", cfg->output_dir, web_files[i]);
        if (copy_file(src_path, dst_path) != 0) return -1;
    }

    /* Copy input WAV as audio.wav */
    snprintf(dst_path, sizeof(dst_path), "%s/audio.wav", cfg->output_dir);
    if (copy_file(cfg->input_path, dst_path) != 0) return -1;

    /* Write data.js */
    snprintf(dst_path, sizeof(dst_path), "%s/data.js", cfg->output_dir);
    if (write_data_js(dst_path, cfg, sample_rate, num_frames,
                      mel_frames, linear_frames, bass_env) != 0)
        return -1;

    return 0;
}

/* -------------------------------------------------------------------
 * Main
 * ---------------------------------------------------------------- */

int main(int argc, char **argv)
{
    config_t cfg;

    MD_set_error_handler(NULL);  /* use default stderr handler */

    if (parse_args(argc, argv, &cfg) != 0)
        return 1;

    /* Read WAV */
    double *samples = NULL;
    double sample_rate = 0.0;
    unsigned num_samples = read_wav(cfg.input_path, &samples, &sample_rate);
    if (num_samples == 0)
        return 1;

    /* Clamp max freq to Nyquist */
    double nyquist = sample_rate / 2.0;
    if (cfg.max_freq > nyquist) {
        fprintf(stderr, "Note: clamping --max-freq to Nyquist (%.0f Hz)\n", nyquist);
        cfg.max_freq = nyquist;
    }

    printf("Input:  %s (%.0f Hz, %u samples, %.1f sec)\n",
           cfg.input_path, sample_rate, num_samples,
           (double)num_samples / sample_rate);

    /* Compute mel frames */
    double *mel_frames = NULL;
    double *bass_env   = NULL;
    unsigned num_frames = compute_mel_frames(samples, num_samples, sample_rate,
                                            &cfg, &mel_frames, &bass_env);
    if (num_frames == 0) {
        free(samples);
        return 1;
    }

    /* Compute linear frames (same frame count, same band count) */
    double *linear_frames = malloc((size_t)num_frames * cfg.num_mels * sizeof(double));
    if (!linear_frames) {
        fprintf(stderr, "Error: allocation failed\n");
        free(samples);
        free(mel_frames);
        free(bass_env);
        return 1;
    }
    compute_linear_frames(samples, num_samples, sample_rate, &cfg,
                          num_frames, linear_frames);
    free(samples);

    printf("Frames: %u (%.0f fps, FFT=%u, %u mel bands)\n",
           num_frames, (double)cfg.fps, cfg.fft_size, cfg.num_mels);

    /* Assemble output */
    if (assemble_output(argv[0], &cfg, sample_rate, num_frames,
                        mel_frames, linear_frames, bass_env) != 0) {
        free(mel_frames);
        free(linear_frames);
        free(bass_env);
        MD_shutdown();
        return 1;
    }

    printf("Output: %s/\n", cfg.output_dir);
    printf("Open in browser: open %s/index.html\n", cfg.output_dir);

    free(mel_frames);
    free(linear_frames);
    free(bass_env);
    MD_shutdown();
    return 0;
}
