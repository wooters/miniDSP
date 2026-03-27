# VAD parameter optimization

## Prerequisites

Install [uv](https://docs.astral.sh/uv/getting-started/installation/).

## Setup

The LibriVAD data is assumed to be here: ~/projects/LibriVAD/Results

```bash
export LIBRIVAD_DIR=~/projects/LibriVAD/Results
```

## Running Optimization

Sample invocation of a small sanity test (narrow conditions, few files):

```bash
uv run python optimize_vad.py librivad ${LIBRIVAD_DIR} \
    --dataset LibriSpeechConcat \
    --split train-clean-100\
    --noises Babble_noise SSN_noise \
    --snrs 10 20 \
    --max-files 10 \
    --n-trials 50 \
    --breakdown
```

Full optimization run:

```bash
uv run python optimize_vad.py librivad ${LIBRIVAD_DIR} \
    --dataset LibriSpeechConcat \
    --split train-clean-100 \
    --n-trials 300 \
    --breakdown \
    --output best_params.json
```
