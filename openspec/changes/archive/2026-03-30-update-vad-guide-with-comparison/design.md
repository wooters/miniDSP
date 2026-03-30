## Context

The VAD guide (`guides/vad.md`) currently covers implementation details, feature formulas, optimization methodology, and API usage. The comparison against ViT-MFCC (small) has been completed and documented in `compare/VAD/README_RESULTS.md`. The guide needs a new section that curates the most useful comparison findings for practitioners deciding between miniDSP VAD and a neural baseline.

The existing guide ends with the "Default parameter optimization" section (including re-tuning guidance) followed by the "API summary". The comparison section should be placed between these — after optimization (which establishes the miniDSP performance ceiling) and before the API reference.

## Goals / Non-Goals

**Goals:**
- Add a self-contained "Comparison with ViT-MFCC baseline" section to `guides/vad.md` that helps practitioners make informed deployment decisions.
- Present headline metrics, per-condition breakdowns, and the dataset-specific optimization experiment in a clear, digestible format.
- Explain what AUC vs F2 reveals about score quality vs threshold tuning — the most non-obvious finding.
- Reference the full comparison tooling for reproducibility.

**Non-Goals:**
- No changes to the comparison code or results themselves.
- No new HTML visualizations or interactive plots for the comparison data.
- No changes to the existing guide sections (features, normalization, state machine, optimization, API).
- Not rewriting or duplicating the full `README_RESULTS.md` — curate and condense for the tutorial audience.

## Decisions

### 1. Place comparison section after optimization, before API summary

**Rationale:** The optimization section establishes the miniDSP F2 ceiling (~0.935). The comparison section naturally follows by showing how that ceiling compares to a neural baseline. The API summary stays last as a reference section.

**Alternative considered:** Placing it as a standalone guide page. Rejected because the comparison is most useful in the context of the VAD tutorial, not as a separate document.

### 2. Condense rather than reproduce README_RESULTS.md

**Rationale:** The guide audience wants actionable takeaways, not raw data dumps. Include the overall results table, a condensed per-SNR table, key noise-type observations (best/worst cases), the dataset-specific optimization headline finding, and the key takeaways. Omit the full per-noise-type table for optimized parameters — reference `compare/VAD/README_RESULTS.md` for the complete breakdown.

### 3. Use existing Markdown table format consistent with the guide

**Rationale:** The guide already uses pipe-delimited Markdown tables (e.g., the state machine table, the optimization results table). The comparison tables should follow the same style for consistency.

## Risks / Trade-offs

- **Stale data risk** → If the comparison is re-run with updated parameters or a different model, the guide numbers will need updating. Mitigated by referencing the results file for canonical numbers and noting the evaluation date.
- **Information density** → Adding ~100–150 lines to an already substantial guide. Mitigated by keeping the section concise and using tables rather than prose for numeric results.
