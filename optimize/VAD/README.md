# VAD parameter optimization

## Prerequisites

Install [uv](https://docs.astral.sh/uv/getting-started/installation/).

## Setup

The LibriVAD project root is assumed to be here: ~/projects/LibriVAD

```bash
export LIBRIVAD_DIR=~/projects/LibriVAD
```

## Running Optimization

Sample invocation of a small sanity test (narrow conditions, few files):

```bash
uv run python optimize_vad.py \
    --n-trials 50 \
    --breakdown \
    librivad ${LIBRIVAD_DIR} \
    --dataset LibriSpeechConcat \
    --split train-clean-100 \
    --noises Babble_noise SSN_noise \
    --snrs 10 20 \
    --max-files 10
```

Full optimization run:

```bash
uv run python optimize_vad.py \
    --n-trials 300 \
    --breakdown \
    --output best_params.json \
    librivad ${LIBRIVAD_DIR} \
    --dataset LibriSpeechConcat \
    --split train-clean-100
```
