## 1. Add comparison section to VAD guide

- [x] 1.1 Add `## Comparison with ViT-MFCC baseline` section to `guides/vad.md` after the "Re-tuning for your data" subsection and before the `## API summary` section, with a horizontal rule separator
- [x] 1.2 Add systems-under-test summary (table or description) covering both systems' type, parameter count, frame rate, training data, and evaluation dataset (LibriVAD test-clean, 702 files)
- [x] 1.3 Add overall results table with F2, precision, recall, AUC (macro), AUC (pooled), and wall time for both systems
- [x] 1.4 Add per-SNR breakdown table showing F2 and AUC for both systems across the six SNR levels, with prose noting the F2 gap consistency (~0.11-0.13) and AUC gap trend (~0.30-0.34)
- [x] 1.5 Add noise-type observations paragraph identifying best cases (babble, SSN) and worst cases (transport, office), with explanation of AUC vs F2 divergence
- [x] 1.6 Add dataset-specific optimization findings subsection: three parameter sets all achieve ~0.935 F2, confirming generalization; F2 ceiling is ~0.935 from tuning alone; remaining gap to ViT requires richer features
- [x] 1.7 Add interpretation guidance explaining AUC vs F2 meaning, threshold sensitivity implications, and practical deployment guidance (embedded/real-time vs compute-available)
- [x] 1.8 Add pointer to `compare/VAD/` for reproduction and `compare/VAD/README_RESULTS.md` for complete breakdowns

## 2. Regenerate llms.txt

- [x] 2.1 Regenerate `llms.txt` and `llms-full.txt` to include the new comparison content
