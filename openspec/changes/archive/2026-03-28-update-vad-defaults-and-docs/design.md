## Context

The VAD module shipped with hand-picked default parameters (equal 0.2 weights, threshold 0.5, onset 3, hangover 15, band 300-3400 Hz). A 300-trial Optuna hyperparameter search against LibriVAD (train-clean-100, all noise types and SNRs) found parameters that improve F2 from 0.837 to 0.933. The optimized values are stored in `optimize/VAD/best_params.json`.

The code change itself is small (one function body), but documentation is the primary deliverable: users and future maintainers need to understand *how* and *why* these defaults were chosen.

## Goals / Non-Goals

**Goals:**
- Replace placeholder defaults with empirically optimized values in `MD_vad_default_params()`.
- Update the `vad-core` spec to be the normative reference for the new defaults.
- Document the optimization methodology, dataset, metric, and results in the tutorial guide, Doxygen header comments, and README.
- Keep all tests passing with the new defaults.

**Non-Goals:**
- Re-running or extending the optimization (the 300-trial run is final for this change).
- Changing the VAD algorithm, API surface, or struct layout.
- Adding runtime parameter loading from JSON — the optimized values are compiled in as constants.
- Optimizing for a different metric (e.g., F1 or equal-error-rate).

## Decisions

### 1. Round parameter values to 6 significant digits

The optimization output has 15+ decimal places. We'll round to 6 significant figures for readability (e.g., `0.723068` not `0.7230676974589997`). The optimization log already shows pre-rounded values in the "C code" snippet — use those.

**Rationale:** 6 digits exceeds the resolution of the optimization (which ran ~2 min/trial with stochastic evaluation). More digits would suggest false precision.

### 2. Use the normalized weights from best_params.json

The raw Optuna weights are unnormalized (sum > 1). `best_params.json` includes `weights_normalized` which sum to ~1.0. Use those, matching the existing convention where default weights sum to 1.0.

**Rationale:** The weighted score feeds into a threshold comparison. Using normalized weights keeps the threshold interpretable and consistent with the existing API contract.

### 3. Document provenance inline in Doxygen and guide

Rather than just stating the new default values, the Doxygen comment for `MD_vad_default_params()` will include a brief note on the optimization source (dataset, metric, trial count, baseline vs. result). The guide will have a dedicated section with full methodology and per-condition breakdown.

**Rationale:** Users tuning VAD for their own domain need to know *what* the defaults were optimized for to decide whether to use them as-is or re-optimize.

### 4. Add optimization methodology section to guides/vad.md

Place a new section "## Default parameter optimization" near the end of the guide (after the algorithm sections, before the API summary). This keeps the tutorial flow intact while making the optimization story discoverable.

### 5. Update README VAD optimization section

The top-level README already has a "VAD optimization" section (added in commit 7bced9d). Update it with the results summary and a pointer to the guide for full methodology.

## Risks / Trade-offs

- **Behavioral change for existing users** — Anyone calling `MD_vad_default_params()` gets different behavior after updating. This is intentional (the old defaults were placeholders) and not an API break, but the documentation should call it out clearly. → *Mitigation:* Note the change prominently in the guide and Doxygen docs.

- **Optimized on one dataset/split** — The parameters were tuned on LibriVAD train-clean-100. Performance on other domains (telephony, far-field, non-English) is unknown. → *Mitigation:* Document the training conditions. The optimization script is already available for users to re-tune on their own data.

- **F2-biased trade-off** — The optimized defaults favor recall (0.98) over precision (0.78). Some applications prefer balanced or precision-biased detection. → *Mitigation:* Document the trade-off explicitly. Users can adjust `threshold` upward for more precision.
