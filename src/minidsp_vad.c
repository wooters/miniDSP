/**
 * @file minidsp_vad.c
 * @brief Voice Activity Detection (VAD) with adaptive feature normalization.
 */

#include "minidsp.h"
#include "minidsp_internal.h"

/* -----------------------------------------------------------------------
 * Feature extraction helpers (static)
 * -----------------------------------------------------------------------*/

/**
 * Spectral entropy: normalize PSD to a probability distribution,
 * return -sum(p * log(p)) / log(num_bins).  Result in [0, 1].
 */
static double compute_spectral_entropy(const double *psd, unsigned num_bins)
{
    double total = 0.0;
    for (unsigned i = 0; i < num_bins; i++)
        total += psd[i];

    if (total <= 0.0)
        return 1.0; /* no energy = maximally unstructured */

    double entropy = 0.0;
    for (unsigned i = 0; i < num_bins; i++) {
        double p = psd[i] / total;
        if (p > 0.0)
            entropy -= p * log(p);
    }

    double max_entropy = log((double)num_bins);
    if (max_entropy <= 0.0)
        return 0.0;

    return entropy / max_entropy;
}

/**
 * Spectral flatness: geometric mean / arithmetic mean of PSD bins.
 * Result in [0, 1].  1.0 = white noise, 0.0 = pure tone.
 */
static double compute_spectral_flatness(const double *psd, unsigned num_bins)
{
    double log_sum = 0.0;
    double arith_sum = 0.0;

    for (unsigned i = 0; i < num_bins; i++) {
        double val = psd[i] > 0.0 ? psd[i] : 1e-30;
        log_sum += log(val);
        arith_sum += psd[i];
    }

    double arith_mean = arith_sum / (double)num_bins;
    if (arith_mean <= 0.0)
        return 1.0; /* no energy = maximally flat */

    double log_geo_mean = log_sum / (double)num_bins;
    double geo_mean = exp(log_geo_mean);

    double flatness = geo_mean / arith_mean;
    if (flatness > 1.0) flatness = 1.0;
    if (flatness < 0.0) flatness = 0.0;

    return flatness;
}

/**
 * Band energy ratio: sum of PSD bins in [band_low_hz, band_high_hz]
 * divided by total PSD sum.  Result in [0, 1].
 */
static double compute_band_energy_ratio(const double *psd, unsigned num_bins,
                                        double sample_rate, unsigned N,
                                        double band_low_hz, double band_high_hz)
{
    double freq_per_bin = sample_rate / (double)N;
    double total = 0.0;
    double band = 0.0;

    for (unsigned i = 0; i < num_bins; i++) {
        double freq = i * freq_per_bin;
        total += psd[i];
        if (freq >= band_low_hz && freq <= band_high_hz)
            band += psd[i];
    }

    if (total <= 0.0)
        return 0.0;

    return band / total;
}

/* -----------------------------------------------------------------------
 * Adaptive normalization helpers (static)
 * -----------------------------------------------------------------------*/

#define RANGE_FLOOR 1e-12

static void update_normalization(MD_vad_state *state, const double *raw)
{
    double alpha = state->params.adaptation_rate;

    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++) {
        if (state->frames_processed == 0) {
            state->feat_min[i] = raw[i];
            state->feat_max[i] = raw[i];
        } else {
            if (raw[i] < state->feat_min[i])
                state->feat_min[i] = state->feat_min[i]
                    + alpha * (raw[i] - state->feat_min[i]);
            if (raw[i] > state->feat_max[i])
                state->feat_max[i] = state->feat_max[i]
                    + alpha * (raw[i] - state->feat_max[i]);
        }
    }
}

static void normalize_features(const MD_vad_state *state,
                               const double *raw, double *norm_out)
{
    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++) {
        double range = state->feat_max[i] - state->feat_min[i];
        if (range < RANGE_FLOOR)
            range = RANGE_FLOOR;

        double val = (raw[i] - state->feat_min[i]) / range;
        if (val < 0.0) val = 0.0;
        if (val > 1.0) val = 1.0;
        norm_out[i] = val;
    }
}

/* -----------------------------------------------------------------------
 * Extract all five raw features from one frame
 * -----------------------------------------------------------------------*/

static void extract_features(const double *signal, unsigned N,
                             double sample_rate,
                             double band_low_hz, double band_high_hz,
                             double *raw_out)
{
    raw_out[MD_VAD_FEAT_ENERGY] = MD_energy(signal, N);
    raw_out[MD_VAD_FEAT_ZCR] = MD_zero_crossing_rate(signal, N);

    unsigned num_bins = N / 2 + 1;
    double psd[num_bins];
    MD_power_spectral_density(signal, N, psd);

    /* Invert entropy and flatness so that higher = more speech-like.
     * Both are naturally high for noise (uniform/flat spectrum) and
     * low for speech (structured harmonics). */
    raw_out[MD_VAD_FEAT_SPECTRAL_ENTROPY] =
        1.0 - compute_spectral_entropy(psd, num_bins);
    raw_out[MD_VAD_FEAT_SPECTRAL_FLATNESS] =
        1.0 - compute_spectral_flatness(psd, num_bins);
    raw_out[MD_VAD_FEAT_BAND_ENERGY_RATIO] =
        compute_band_energy_ratio(psd, num_bins, sample_rate, N,
                                  band_low_hz, band_high_hz);
}

/* -----------------------------------------------------------------------
 * Public API
 * -----------------------------------------------------------------------*/

void MD_vad_default_params(MD_vad_params *params)
{
    MD_CHECK_VOID(params != NULL, MD_ERR_NULL_POINTER, "params is NULL");

    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++)
        params->weights[i] = 0.2;

    params->threshold       = 0.5;
    params->onset_frames    = 3;
    params->hangover_frames = 15;
    params->adaptation_rate = 0.01;
    params->band_low_hz     = 300.0;
    params->band_high_hz    = 3400.0;
}

void MD_vad_init(MD_vad_state *state, const MD_vad_params *params)
{
    MD_CHECK_VOID(state != NULL, MD_ERR_NULL_POINTER, "state is NULL");

    if (params != NULL) {
        state->params = *params;
    } else {
        MD_vad_default_params(&state->params);
    }

    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++) {
        state->feat_min[i] = 1e30;
        state->feat_max[i] = -1e30;
    }

    state->onset_counter    = 0;
    state->hangover_counter = 0;
    state->current_decision = 0;
    state->frames_processed = 0;
}

void MD_vad_calibrate(MD_vad_state *state, const double *signal,
                      unsigned N, double sample_rate)
{
    MD_CHECK_VOID(state != NULL, MD_ERR_NULL_POINTER, "state is NULL");
    MD_CHECK_VOID(signal != NULL, MD_ERR_NULL_POINTER, "signal is NULL");
    MD_CHECK_VOID(N >= 2, MD_ERR_INVALID_SIZE, "N must be >= 2");

    double raw[MD_VAD_NUM_FEATURES];
    extract_features(signal, N, sample_rate,
                     state->params.band_low_hz, state->params.band_high_hz,
                     raw);

    update_normalization(state, raw);
    state->frames_processed++;
}

int MD_vad_process_frame(MD_vad_state *state, const double *signal,
                         unsigned N, double sample_rate,
                         double *score_out, double *features_out)
{
    MD_CHECK(state != NULL, MD_ERR_NULL_POINTER, "state is NULL", 0);
    MD_CHECK(signal != NULL, MD_ERR_NULL_POINTER, "signal is NULL", 0);
    MD_CHECK(N >= 2, MD_ERR_INVALID_SIZE, "N must be >= 2", 0);

    /* 1. Extract raw features */
    double raw[MD_VAD_NUM_FEATURES];
    extract_features(signal, N, sample_rate,
                     state->params.band_low_hz, state->params.band_high_hz,
                     raw);

    /* 2. Update adaptive normalization */
    update_normalization(state, raw);

    /* 3. Normalize features to [0, 1] */
    double norm[MD_VAD_NUM_FEATURES];
    normalize_features(state, raw, norm);

    /* 4. Compute weighted score */
    double score = 0.0;
    for (int i = 0; i < MD_VAD_NUM_FEATURES; i++)
        score += state->params.weights[i] * norm[i];

    /* 5. Apply state machine */
    if (score >= state->params.threshold) {
        state->onset_counter++;
        if (state->current_decision == 0) {
            if (state->onset_counter >= state->params.onset_frames) {
                state->current_decision = 1;
                state->hangover_counter = state->params.hangover_frames;
            }
        } else {
            state->hangover_counter = state->params.hangover_frames;
        }
    } else {
        if (state->current_decision == 1) {
            if (state->hangover_counter > 0) {
                state->hangover_counter--;
            }
            if (state->hangover_counter == 0) {
                state->current_decision = 0;
                state->onset_counter = 0;
            }
        } else {
            state->onset_counter = 0;
        }
    }

    /* 6. Write optional outputs */
    if (score_out != NULL)
        *score_out = score;
    if (features_out != NULL) {
        for (int i = 0; i < MD_VAD_NUM_FEATURES; i++)
            features_out[i] = norm[i];
    }

    /* 7. Increment counter, return decision */
    state->frames_processed++;
    return state->current_decision;
}
