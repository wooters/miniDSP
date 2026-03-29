#!/usr/bin/env python3
"""Hyperparameter optimization for miniDSP VAD using Optuna + pyminidsp.

Works directly with the LibriVAD dataset structure:
  - Audio: {librivad_root}/Results/{dataset}/{split}/{noise}/{snr}/{file}.wav
  - Labels: {librivad_root}/Files/Labels/{dataset}/{split}/{spk}/{ch}/{file}.npy

Labels are sample-level binary arrays (16 kHz, 1=speech, 0=silence) from
LibriSpeech forced alignments.  They are downsampled to frame-level via
majority voting to match the VAD's frame rate.

Also supports plain Audacity label files for hand-labeled audio.

Usage (LibriVAD):
    python optimize_vad.py librivad /path/to/LibriVAD \
        --noises Babble_noise SSN_noise \
        --snrs 10 20 \
        --dataset LibriSpeechConcat \
        --split test-clean \
        --n-trials 300

Usage (Audacity labels):
    python optimize_vad.py audacity \
        --audio a.wav b.wav \
        --labels a_labels.txt b_labels.txt \
        --n-trials 300

Dependencies:
    pip install pyminidsp optuna soundfile numpy scipy
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import numpy.typing as npt
import optuna
import soundfile as sf

from pyminidsp import VAD


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EvalItem:
    """One audio file paired with its frame-level ground truth."""

    audio: npt.NDArray[np.float64]
    sample_rate: float
    targets: npt.NDArray[np.int32]
    name: str
    noise: str = ""
    snr: int = 0


@dataclass(frozen=True)
class Segment:
    """A labeled time segment (for Audacity format)."""

    start: float
    end: float
    label: str


# ---------------------------------------------------------------------------
# LibriVAD loading
# ---------------------------------------------------------------------------

def resolve_label_path(
    wav_path: Path,
    labels_root: Path,
    dataset: str,
    split: str,
) -> Path:
    """Map a LibriVAD noisy WAV path to its sample-level label .npy file.

    Filename like ``1089-134686-0000.wav`` maps to
    ``{labels_root}/{dataset}/{split}/1089/134686/1089-134686-0000.npy``.

    For LibriSpeechConcat files (containing ``_+_``), the label lives under
    the LibriSpeechConcat label tree with the full concatenated filename.

    Args:
        wav_path: Path to the noisy WAV file.
        labels_root: Root of the Labels directory tree.
        dataset: ``"LibriSpeech"`` or ``"LibriSpeechConcat"``.
        split: e.g. ``"test-clean"``, ``"dev-clean"``, ``"train-clean-100"``.

    Returns:
        Path to the corresponding ``.npy`` label file.
    """
    stem = wav_path.stem
    parts = stem.split("_+_")[0].split("-")  # always use first part for spk/ch
    speaker = parts[0]
    chapter = parts[1]
    return labels_root / dataset / split / speaker / chapter / f"{stem}.npy"


def sample_labels_to_frame_labels(
    sample_labels: npt.NDArray[np.int16],
    frame_len_samples: int,
) -> npt.NDArray[np.int32]:
    """Downsample sample-level binary labels to frame-level via majority vote.

    Args:
        sample_labels: 1-D array of 0/1 at sample rate (e.g. 16 kHz).
        frame_len_samples: Number of samples per VAD frame.

    Returns:
        int32 array of frame-level labels (1=speech, 0=nonspeech).
    """
    num_frames = len(sample_labels) // frame_len_samples
    if num_frames == 0:
        return np.array([], dtype=np.int32)

    trimmed = sample_labels[: num_frames * frame_len_samples]
    frames = trimmed.reshape(num_frames, frame_len_samples)
    return (frames.mean(axis=1) >= 0.5).astype(np.int32)


def discover_librivad_files(
    results_root: Path,
    dataset: str,
    split: str,
    noises: list[str] | None,
    snrs: list[int] | None,
    max_files_per_condition: int | None = None,
) -> list[tuple[Path, str, int]]:
    """Find noisy WAV files in the LibriVAD Results tree.

    Args:
        results_root: Path to the ``Results/`` directory.
        dataset: ``"LibriSpeech"`` or ``"LibriSpeechConcat"``.
        split: e.g. ``"test-clean"``.
        noises: Noise types to include, or None for all.
        snrs: SNR levels to include, or None for all.
        max_files_per_condition: Cap files per noise/SNR pair (for speed).

    Returns:
        List of (wav_path, noise_name, snr_value) tuples.
    """
    dataset_dir = results_root / dataset / split
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset directory not found: {dataset_dir}")

    noise_dirs = (
        sorted(d for d in dataset_dir.iterdir() if d.is_dir())
        if noises is None
        else [dataset_dir / n for n in noises]
    )

    files: list[tuple[Path, str, int]] = []
    for noise_dir in noise_dirs:
        if not noise_dir.is_dir():
            continue
        noise_name = noise_dir.name

        snr_dirs = (
            sorted(d for d in noise_dir.iterdir() if d.is_dir())
            if snrs is None
            else [noise_dir / str(s) for s in snrs]
        )

        for snr_dir in snr_dirs:
            if not snr_dir.is_dir():
                continue
            snr_val = int(snr_dir.name)

            wavs = sorted(snr_dir.glob("*.wav"))
            if max_files_per_condition is not None:
                wavs = wavs[:max_files_per_condition]

            for w in wavs:
                files.append((w, noise_name, snr_val))

    return files


def load_librivad_dataset(
    librivad_root: Path,
    dataset: str,
    split: str,
    noises: list[str] | None,
    snrs: list[int] | None,
    frame_len_ms: float = 20.0,
    max_files_per_condition: int | None = None,
) -> list[EvalItem]:
    """Load a LibriVAD evaluation dataset.

    Args:
        librivad_root: Root of the LibriVAD project (contains Results/ and Files/).
        dataset: ``"LibriSpeech"`` or ``"LibriSpeechConcat"``.
        split: e.g. ``"test-clean"``.
        noises: Noise types to include, or None for all.
        snrs: SNR levels to include, or None for all.
        frame_len_ms: VAD frame length in ms.
        max_files_per_condition: Cap files per noise/SNR pair.

    Returns:
        List of EvalItem objects.
    """
    results_root = librivad_root / "Results"
    labels_root = librivad_root / "Files" / "Labels"

    discovered = discover_librivad_files(
        results_root, dataset, split, noises, snrs, max_files_per_condition,
    )
    if not discovered:
        raise RuntimeError(
            f"No WAV files found under {results_root / dataset / split}. "
            f"Check --dataset, --split, --noises, --snrs."
        )

    items: list[EvalItem] = []
    skipped = 0

    for wav_path, noise_name, snr_val in discovered:
        label_path = resolve_label_path(wav_path, labels_root, dataset, split)
        if not label_path.exists():
            skipped += 1
            continue

        audio, sr = sf.read(str(wav_path), dtype="float64")
        if audio.ndim > 1:
            audio = audio.mean(axis=1)

        sample_labels = np.load(label_path).astype(np.int16)

        frame_len_samples = int(sr * frame_len_ms / 1000.0)

        # Labels and audio may differ slightly in length; truncate to shorter.
        min_samples = min(len(audio), len(sample_labels))
        audio = audio[:min_samples]
        sample_labels = sample_labels[:min_samples]

        targets = sample_labels_to_frame_labels(sample_labels, frame_len_samples)
        num_frames = len(targets)

        audio = audio[: num_frames * frame_len_samples]

        items.append(EvalItem(
            audio=audio,
            sample_rate=sr,
            targets=targets,
            name=wav_path.stem,
            noise=noise_name,
            snr=snr_val,
        ))

    if skipped > 0:
        print(f"  Warning: skipped {skipped} files with missing labels.")

    return items


# ---------------------------------------------------------------------------
# Audacity label loading (kept for non-LibriVAD use)
# ---------------------------------------------------------------------------

def load_audacity_labels(path: Path) -> list[Segment]:
    """Load an Audacity label file (tab-separated: start, end, label).

    Args:
        path: Path to the label file.

    Returns:
        List of Segment objects sorted by start time.
    """
    segments: list[Segment] = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            segments.append(Segment(
                start=float(parts[0]),
                end=float(parts[1]),
                label=parts[2].strip().lower(),
            ))
    return sorted(segments, key=lambda s: s.start)


def audacity_labels_to_frame_targets(
    segments: list[Segment],
    num_frames: int,
    frame_dur: float,
) -> npt.NDArray[np.int32]:
    """Convert time-based labels to per-frame binary targets.

    A frame is speech (1) if its midpoint falls within a speech segment.

    Args:
        segments: Labeled segments.
        num_frames: Total frames.
        frame_dur: Duration of each frame in seconds.

    Returns:
        int32 array of per-frame labels.
    """
    targets = np.zeros(num_frames, dtype=np.int32)
    for i in range(num_frames):
        frame_mid = i * frame_dur + frame_dur / 2.0
        for seg in segments:
            if seg.start > frame_mid:
                break
            if seg.label == "speech" and seg.start <= frame_mid < seg.end:
                targets[i] = 1
                break
    return targets


def load_audacity_dataset(
    audio_paths: list[Path],
    label_paths: list[Path],
    frame_len_ms: float = 20.0,
) -> list[EvalItem]:
    """Load audio/label pairs using Audacity label files.

    Args:
        audio_paths: Paths to WAV files.
        label_paths: Corresponding Audacity label files.
        frame_len_ms: Frame length in ms.

    Returns:
        List of EvalItem objects.
    """
    items: list[EvalItem] = []
    for audio_path, label_path in zip(audio_paths, label_paths):
        audio, sr = sf.read(str(audio_path), dtype="float64")
        if audio.ndim > 1:
            audio = audio.mean(axis=1)

        frame_len_samples = int(sr * frame_len_ms / 1000.0)
        num_frames = len(audio) // frame_len_samples
        frame_dur = frame_len_samples / sr

        segments = load_audacity_labels(label_path)
        targets = audacity_labels_to_frame_targets(segments, num_frames, frame_dur)

        items.append(EvalItem(
            audio=audio,
            sample_rate=sr,
            targets=targets,
            name=audio_path.stem,
        ))
    return items


# ---------------------------------------------------------------------------
# Evaluation
# ---------------------------------------------------------------------------

def compute_fbeta(
    predictions: npt.NDArray[np.int32],
    targets: npt.NDArray[np.int32],
    beta: float = 2.0,
) -> dict[str, float]:
    """Compute F-beta score and supporting metrics.

    Args:
        predictions: Binary predictions (0 or 1).
        targets: Binary ground truth (0 or 1).
        beta: F-beta parameter. beta > 1 favors recall.

    Returns:
        Dict with f_beta, precision, recall, tp, fp, fn, tn.
    """
    tp = int(np.sum((predictions == 1) & (targets == 1)))
    fp = int(np.sum((predictions == 1) & (targets == 0)))
    fn = int(np.sum((predictions == 0) & (targets == 1)))
    tn = int(np.sum((predictions == 0) & (targets == 0)))

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0

    beta_sq = beta * beta
    denom = beta_sq * precision + recall
    f_beta = (1.0 + beta_sq) * precision * recall / denom if denom > 0.0 else 0.0

    return {
        "f_beta": f_beta,
        "precision": precision,
        "recall": recall,
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
    }


def _eval_one(args: tuple) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.int32]]:
    """Worker: run VAD on a single audio file and return (predictions, targets).

    Args:
        args: Tuple of (audio, sample_rate, targets, vad_kwargs, frame_len_ms).

    Returns:
        (predictions, targets) arrays.
    """
    audio, sample_rate, targets, vad_kwargs, frame_len_ms = args
    frame_len_samples = int(sample_rate * frame_len_ms / 1000.0)
    vad = VAD(**vad_kwargs)
    decisions, _, _ = vad.process(audio, sample_rate, frame_len_samples)
    return decisions, targets


def evaluate_params(
    dataset: list[EvalItem],
    *,
    weights: list[float],
    threshold: float,
    onset_frames: int,
    hangover_frames: int,
    adaptation_rate: float,
    band_low_hz: float,
    band_high_hz: float,
    frame_len_ms: float = 20.0,
    beta: float = 2.0,
    executor: ProcessPoolExecutor | None = None,
) -> dict[str, float]:
    """Run VAD with given params on all items, return aggregate F-beta.

    Args:
        dataset: List of EvalItem objects.
        weights: Feature weights (length 5).
        threshold: Decision threshold.
        onset_frames: Onset frame count.
        hangover_frames: Hangover frame count.
        adaptation_rate: EMA adaptation rate.
        band_low_hz: Speech band lower bound (Hz).
        band_high_hz: Speech band upper bound (Hz).
        frame_len_ms: Frame length in ms.
        beta: F-beta parameter.
        executor: Process pool for parallel evaluation. Serial if None.

    Returns:
        Aggregate metrics dict.
    """
    vad_kwargs = dict(
        weights=weights,
        threshold=threshold,
        onset_frames=onset_frames,
        hangover_frames=hangover_frames,
        adaptation_rate=adaptation_rate,
        band_low_hz=band_low_hz,
        band_high_hz=band_high_hz,
    )

    work = [
        (item.audio, item.sample_rate, item.targets, vad_kwargs, frame_len_ms)
        for item in dataset
    ]

    if executor is not None:
        results = list(executor.map(_eval_one, work))
    else:
        results = [_eval_one(w) for w in work]

    preds = np.concatenate([r[0] for r in results])
    targets = np.concatenate([r[1] for r in results])
    return compute_fbeta(preds, targets, beta=beta)


def evaluate_per_condition(
    dataset: list[EvalItem],
    *,
    weights: list[float],
    threshold: float,
    onset_frames: int,
    hangover_frames: int,
    adaptation_rate: float,
    band_low_hz: float,
    band_high_hz: float,
    frame_len_ms: float = 20.0,
    beta: float = 2.0,
    executor: ProcessPoolExecutor | None = None,
) -> dict[str, dict[str, float]]:
    """Evaluate per noise type and SNR, returning a breakdown.

    Args:
        dataset: List of EvalItem objects (with noise/snr fields set).
        weights: Feature weights (length 5).
        threshold: Decision threshold.
        onset_frames: Onset frame count.
        hangover_frames: Hangover frame count.
        adaptation_rate: EMA adaptation rate.
        band_low_hz: Speech band lower bound (Hz).
        band_high_hz: Speech band upper bound (Hz).
        frame_len_ms: Frame length in ms.
        beta: F-beta parameter.
        executor: Process pool for parallel evaluation. Serial if None.

    Returns:
        Dict mapping ``"{noise}@{snr}dB"`` to metrics dicts.
    """
    buckets: dict[str, list[EvalItem]] = defaultdict(list)
    for item in dataset:
        key = f"{item.noise}@{item.snr}dB"
        buckets[key].append(item)

    results: dict[str, dict[str, float]] = {}
    for key, items in sorted(buckets.items()):
        results[key] = evaluate_params(
            items,
            weights=weights,
            threshold=threshold,
            onset_frames=onset_frames,
            hangover_frames=hangover_frames,
            adaptation_rate=adaptation_rate,
            band_low_hz=band_low_hz,
            band_high_hz=band_high_hz,
            frame_len_ms=frame_len_ms,
            beta=beta,
            executor=executor,
        )
    return results


# ---------------------------------------------------------------------------
# Optuna objective
# ---------------------------------------------------------------------------

def make_objective(
    dataset: list[EvalItem],
    frame_len_ms: float = 20.0,
    beta: float = 2.0,
    executor: ProcessPoolExecutor | None = None,
) -> callable:
    """Create an Optuna objective closed over the dataset.

    Args:
        dataset: Evaluation dataset.
        frame_len_ms: Frame length in ms.
        beta: F-beta parameter.
        executor: Process pool for parallel evaluation. Serial if None.

    Returns:
        Objective function for Optuna.
    """

    def objective(trial: optuna.Trial) -> float:
        raw_weights = [
            trial.suggest_float(f"w_{name}", 0.01, 5.0, log=True)
            for name in ["energy", "zcr", "entropy", "flatness", "band_ratio"]
        ]
        total = sum(raw_weights)
        weights = [w / total for w in raw_weights]

        threshold = trial.suggest_float("threshold", 0.15, 0.85)
        onset_frames = trial.suggest_int("onset_frames", 1, 10)
        hangover_frames = trial.suggest_int("hangover_frames", 3, 40)
        adaptation_rate = trial.suggest_float("adaptation_rate", 0.001, 0.1, log=True)
        band_low_hz = trial.suggest_float("band_low_hz", 80.0, 500.0)
        band_high_hz = trial.suggest_float("band_high_hz", 2000.0, 5000.0)

        metrics = evaluate_params(
            dataset,
            weights=weights,
            threshold=threshold,
            onset_frames=onset_frames,
            hangover_frames=hangover_frames,
            adaptation_rate=adaptation_rate,
            band_low_hz=band_low_hz,
            band_high_hz=band_high_hz,
            frame_len_ms=frame_len_ms,
            beta=beta,
            executor=executor,
        )

        trial.set_user_attr("precision", metrics["precision"])
        trial.set_user_attr("recall", metrics["recall"])
        trial.set_user_attr("weights_normalized", weights)

        return metrics["f_beta"]

    return objective


# ---------------------------------------------------------------------------
# Results formatting
# ---------------------------------------------------------------------------

def format_c_defaults(trial: optuna.trial.FrozenTrial) -> str:
    """Format best trial as C code for ``MD_vad_default_params``.

    Args:
        trial: The best Optuna trial.

    Returns:
        C code string.
    """
    weights = trial.user_attrs["weights_normalized"]
    p = trial.params

    lines = [
        "/* Optimized VAD parameters (F2-optimized, recall-biased) */",
        f"params->weights[0] = {weights[0]:.6f};  /* energy */",
        f"params->weights[1] = {weights[1]:.6f};  /* zcr */",
        f"params->weights[2] = {weights[2]:.6f};  /* spectral_entropy */",
        f"params->weights[3] = {weights[3]:.6f};  /* spectral_flatness */",
        f"params->weights[4] = {weights[4]:.6f};  /* band_energy_ratio */",
        f"params->threshold       = {p['threshold']:.6f};",
        f"params->onset_frames    = {p['onset_frames']};",
        f"params->hangover_frames = {p['hangover_frames']};",
        f"params->adaptation_rate = {p['adaptation_rate']:.6f};",
        f"params->band_low_hz     = {p['band_low_hz']:.1f};",
        f"params->band_high_hz    = {p['band_high_hz']:.1f};",
    ]
    return "\n".join(lines)


def format_python_constructor(trial: optuna.trial.FrozenTrial) -> str:
    """Format best trial as a pyminidsp ``VAD()`` call.

    Args:
        trial: The best Optuna trial.

    Returns:
        Python code string.
    """
    weights = trial.user_attrs["weights_normalized"]
    p = trial.params

    w_str = ", ".join(f"{w:.6f}" for w in weights)
    return (
        f"VAD(\n"
        f"    weights=[{w_str}],\n"
        f"    threshold={p['threshold']:.6f},\n"
        f"    onset_frames={p['onset_frames']},\n"
        f"    hangover_frames={p['hangover_frames']},\n"
        f"    adaptation_rate={p['adaptation_rate']:.6f},\n"
        f"    band_low_hz={p['band_low_hz']:.1f},\n"
        f"    band_high_hz={p['band_high_hz']:.1f},\n"
        f")"
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser with subcommands."""
    parser = argparse.ArgumentParser(
        description="Optimize miniDSP VAD hyperparameters.",
    )
    parser.add_argument(
        "--n-trials", type=int, default=300,
        help="Number of Optuna trials (default: 300)",
    )
    parser.add_argument(
        "--frame-len-ms", type=float, default=20.0,
        help="Frame length in ms (default: 20.0)",
    )
    parser.add_argument(
        "--beta", type=float, default=2.0,
        help="F-beta parameter; >1 favors recall (default: 2.0)",
    )
    parser.add_argument(
        "--output", type=Path, default=None,
        help="Save best params as JSON",
    )
    parser.add_argument(
        "--breakdown", action="store_true",
        help="Print per-condition (noise/SNR) results for best params",
    )
    parser.add_argument(
        "--workers", type=int, default=None,
        help="Number of parallel workers for per-file evaluation "
             "(default: all CPU cores, 0 = serial)",
    )

    sub = parser.add_subparsers(dest="mode", required=True)

    # -- LibriVAD subcommand --
    lv = sub.add_parser("librivad", help="Use LibriVAD dataset")
    lv.add_argument("root", type=Path, help="LibriVAD project root directory")
    lv.add_argument(
        "--dataset", default="LibriSpeechConcat",
        choices=["LibriSpeech", "LibriSpeechConcat"],
        help="Which dataset variant (default: LibriSpeechConcat)",
    )
    lv.add_argument(
        "--split", default="test-clean",
        help="Which split (default: test-clean)",
    )
    lv.add_argument(
        "--noises", nargs="+", default=None,
        help="Noise types to include (default: all)",
    )
    lv.add_argument(
        "--snrs", nargs="+", type=int, default=None,
        help="SNR levels to include (default: all)",
    )
    lv.add_argument(
        "--max-files", type=int, default=None,
        help="Max files per noise/SNR condition (for faster iteration)",
    )

    # -- Audacity subcommand --
    au = sub.add_parser("audacity", help="Use hand-labeled Audacity files")
    au.add_argument("--audio", nargs="+", type=Path, required=True)
    au.add_argument("--labels", nargs="+", type=Path, required=True)

    return parser


def run_baseline(
    dataset: list[EvalItem],
    frame_len_ms: float,
    beta: float,
    executor: ProcessPoolExecutor | None = None,
) -> dict[str, float]:
    """Evaluate the current C library defaults as a baseline.

    Args:
        dataset: Evaluation dataset.
        frame_len_ms: Frame length in ms.
        beta: F-beta parameter.
        executor: Process pool for parallel evaluation. Serial if None.

    Returns:
        Metrics dict for the default params.
    """
    return evaluate_params(
        dataset,
        weights=[0.2, 0.2, 0.2, 0.2, 0.2],
        threshold=0.5,
        onset_frames=3,
        hangover_frames=15,
        adaptation_rate=0.01,
        band_low_hz=300.0,
        band_high_hz=3400.0,
        frame_len_ms=frame_len_ms,
        beta=beta,
        executor=executor,
    )


def main(argv: list[str] | None = None) -> None:
    """Entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)

    # --- Load dataset ---
    if args.mode == "librivad":
        print(f"Loading LibriVAD: {args.root}")
        print(f"  dataset={args.dataset}  split={args.split}")
        print(f"  noises={args.noises or 'all'}  snrs={args.snrs or 'all'}")
        dataset = load_librivad_dataset(
            librivad_root=args.root,
            dataset=args.dataset,
            split=args.split,
            noises=args.noises,
            snrs=args.snrs,
            frame_len_ms=args.frame_len_ms,
            max_files_per_condition=args.max_files,
        )
    else:
        if len(args.audio) != len(args.labels):
            parser.error("Must provide same number of --audio and --labels files.")
        print(f"Loading {len(args.audio)} audio/label pairs...")
        dataset = load_audacity_dataset(args.audio, args.labels, args.frame_len_ms)

    # --- Summary ---
    total_frames = sum(len(item.targets) for item in dataset)
    speech_frames = sum(int(item.targets.sum()) for item in dataset)
    total_secs = total_frames * args.frame_len_ms / 1000.0
    print(f"\nDataset: {len(dataset)} files, {total_frames} frames "
          f"({total_secs:.1f}s), {speech_frames} speech "
          f"({100*speech_frames/total_frames:.1f}%)")

    # --- Set up parallel executor ---
    n_workers = args.workers
    if n_workers is None:
        n_workers = os.cpu_count() or 1
    executor = ProcessPoolExecutor(max_workers=n_workers) if n_workers > 0 else None
    if executor is not None:
        print(f"\nParallel evaluation: {n_workers} workers")

    # --- Baseline ---
    default_metrics = run_baseline(dataset, args.frame_len_ms, args.beta,
                                   executor=executor)
    print(f"\nBaseline (current defaults):")
    print(f"  F{args.beta:.0f}={default_metrics['f_beta']:.4f}  "
          f"P={default_metrics['precision']:.4f}  "
          f"R={default_metrics['recall']:.4f}")

    # --- Optimize ---
    print(f"\nRunning {args.n_trials} trials (F{args.beta:.0f} optimization)...")
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="maximize",
        sampler=optuna.samplers.TPESampler(seed=42),
    )
    objective = make_objective(dataset, args.frame_len_ms, args.beta,
                               executor=executor)
    study.optimize(objective, n_trials=args.n_trials, show_progress_bar=True)

    # --- Report ---
    best = study.best_trial
    print(f"\n{'='*60}")
    print(f"Best trial #{best.number}:")
    print(f"  F{args.beta:.0f}={best.value:.4f}  "
          f"P={best.user_attrs['precision']:.4f}  "
          f"R={best.user_attrs['recall']:.4f}")

    improvement = best.value - default_metrics["f_beta"]
    print(f"  Improvement over defaults: {improvement:+.4f} F{args.beta:.0f}")

    print(f"\n--- C code for MD_vad_default_params ---")
    print(format_c_defaults(best))
    print(f"\n--- Python constructor ---")
    print(format_python_constructor(best))

    # --- Per-condition breakdown ---
    if args.breakdown and args.mode == "librivad":
        weights = best.user_attrs["weights_normalized"]
        p = best.params
        breakdown = evaluate_per_condition(
            dataset,
            weights=weights,
            threshold=p["threshold"],
            onset_frames=p["onset_frames"],
            hangover_frames=p["hangover_frames"],
            adaptation_rate=p["adaptation_rate"],
            band_low_hz=p["band_low_hz"],
            band_high_hz=p["band_high_hz"],
            frame_len_ms=args.frame_len_ms,
            beta=args.beta,
            executor=executor,
        )
        print(f"\n--- Per-condition breakdown (best params) ---")
        print(f"{'Condition':<30s} {'F'+str(int(args.beta)):>6s} {'Prec':>6s} {'Rec':>6s}")
        print("-" * 50)
        for cond, m in breakdown.items():
            print(f"{cond:<30s} {m['f_beta']:6.4f} {m['precision']:6.4f} {m['recall']:6.4f}")

    # --- Save JSON ---
    if args.output:
        result = {
            "best_f_beta": best.value,
            "beta": args.beta,
            "precision": best.user_attrs["precision"],
            "recall": best.user_attrs["recall"],
            "params": best.params,
            "weights_normalized": best.user_attrs["weights_normalized"],
            "baseline_f_beta": default_metrics["f_beta"],
            "n_trials": args.n_trials,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"\nResults saved to {args.output}")

    if executor is not None:
        executor.shutdown()


if __name__ == "__main__":
    main()
