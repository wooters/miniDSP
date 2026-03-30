#!/usr/bin/env python3
"""Compare miniDSP VAD against ViT-MFCC baseline on LibriVAD.

Evaluates both systems on the same test data and reports F-beta, precision,
and recall side by side.

Both systems operate at their native frame rates:
  - miniDSP VAD:  20 ms frames (50 fps)
  - ViT-MFCC:     10 ms frames (100 fps)

Ground-truth labels are downsampled from sample-level to each system's frame
rate independently via majority voting.

Usage:
    uv run python compare_vad.py \\
        --librivad-root ~/projects/LibriVAD \\
        --breakdown

Dependencies:
    See pyproject.toml (install via ``uv sync``).
    Also requires the LibriVAD repo with Eval-ViT-MFCC/ and prepared test data.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import numpy.typing as npt
import soundfile as sf
import torch
from sklearn.metrics import roc_auc_score
from torch import nn, einsum
from einops import rearrange
from einops.layers.torch import Rearrange

from pyminidsp import VAD


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HF_REPO_ID = "LibriVAD/LibriVAD"
HF_REPO_TYPE = "dataset"
HF_MODEL_PATH = (
    "Eval-ViT-MFCC/small-models/LibriSpeechConcat_vit_MFCC/model_epoch_50.model"
)
CACHE_DIR = Path(__file__).resolve().parent / ".model_cache"

# ViT-MFCC model hyperparameters (fixed across all model "sizes")
VIT_INPUT_RES = [39, 100]
VIT_PATCH_RES = [39, 1]
VIT_NUM_CLASSES = 2
VIT_DIM = 192
VIT_DEPTH = 12
VIT_HEADS = 3
VIT_MLP_DIM = 768
VIT_CHANNELS = 1
VIT_DIM_HEAD = VIT_DIM // VIT_HEADS  # 64
VIT_DROPOUT = 0.1
VIT_EMB_DROPOUT = 0.1
VIT_PRE_NORM = "False"
VIT_SEQUENCE_LENGTH = 100

# MFCC extraction parameters (must match ViT training)
MFCC_WINLEN = 0.025   # 25 ms
MFCC_OVRLEN = 0.01    # 10 ms shift
MFCC_PRE_COEF = 0.97
MFCC_NFILTER = 24
MFCC_NFFT = 512
MFCC_NO_COEF = 12

# miniDSP frame rate
MINIDSP_FRAME_MS = 20.0


# ---------------------------------------------------------------------------
# KWT (Keyword Transformer) model — inlined from LibriVAD KWT3_vad.py
# because KWT3_vad.py has top-level script code that prevents clean import.
# Original: https://github.com/wdjose/keyword-transformer
# Copyright (c) 2021 Williard Joshua Jose, MIT License
# ---------------------------------------------------------------------------

class _PreNorm(nn.Module):
    def __init__(self, dim, fn):
        super().__init__()
        self.norm = nn.LayerNorm(dim)
        self.fn = fn

    def forward(self, x, **kwargs):
        return self.fn(self.norm(x), **kwargs)


class _PostNorm(nn.Module):
    def __init__(self, dim, fn):
        super().__init__()
        self.norm = nn.LayerNorm(dim)
        self.fn = fn

    def forward(self, x, **kwargs):
        return self.norm(self.fn(x, **kwargs))


class _FeedForward(nn.Module):
    def __init__(self, dim, hidden_dim, dropout=0.0):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, dim),
            nn.Dropout(dropout),
        )

    def forward(self, x):
        return self.net(x)


class _Attention(nn.Module):
    def __init__(self, dim, heads=8, dim_head=64, dropout=0.0):
        super().__init__()
        inner_dim = dim_head * heads
        project_out = not (heads == 1 and dim_head == dim)
        self.heads = heads
        self.scale = dim_head ** -0.5
        self.attend = nn.Softmax(dim=-1)
        self.to_qkv = nn.Linear(dim, inner_dim * 3, bias=False)
        self.to_out = (
            nn.Sequential(nn.Linear(inner_dim, dim), nn.Dropout(dropout))
            if project_out
            else nn.Identity()
        )

    def forward(self, x):
        b, n, _, h = *x.shape, self.heads
        qkv = self.to_qkv(x).chunk(3, dim=-1)
        q, k, v = map(lambda t: rearrange(t, "b n (h d) -> b h n d", h=h), qkv)
        dots = einsum("b h i d, b h j d -> b h i j", q, k) * self.scale
        attn = self.attend(dots)
        out = einsum("b h i j, b h j d -> b h i d", attn, v)
        out = rearrange(out, "b h n d -> b n (h d)")
        return self.to_out(out)


class _Transformer(nn.Module):
    def __init__(self, dim, depth, heads, dim_head, mlp_dim, pre_norm=True, dropout=0.0):
        super().__init__()
        self.layers = nn.ModuleList([])
        P_Norm = _PreNorm if pre_norm else _PostNorm
        for _ in range(depth):
            self.layers.append(nn.ModuleList([
                P_Norm(dim, _Attention(dim, heads=heads, dim_head=dim_head, dropout=dropout)),
                P_Norm(dim, _FeedForward(dim, mlp_dim, dropout=dropout)),
            ]))

    def forward(self, x):
        for attn, ff in self.layers:
            x = attn(x) + x
            x = ff(x) + x
        return x


class KWT(nn.Module):
    """Keyword Transformer for frame-level VAD classification."""

    def __init__(self, input_res, patch_res, num_classes, dim, depth, heads,
                 mlp_dim, channels, dim_head, dropout, emb_dropout, pre_norm=True, **kwargs):
        super().__init__()
        self.num_patches = int(input_res[0] / patch_res[0] * input_res[1] / patch_res[1])
        self.patch_dim = channels * patch_res[0] * patch_res[1]
        self.to_patch_embedding = nn.Sequential(
            Rearrange("b c (h p1) (w p2) -> b (h w) (p1 p2 c)", p1=patch_res[0], p2=patch_res[1]),
            nn.Linear(self.patch_dim, dim),
        )
        self.pos_embedding = nn.Parameter(torch.randn(1, self.num_patches, dim))
        self.dropout = nn.Dropout(emb_dropout)
        self.transformer = _Transformer(dim, depth, heads, dim_head, mlp_dim, pre_norm, dropout)
        self.to_latent = nn.Identity()
        self.mlp_head = nn.Sequential(nn.LayerNorm(dim), nn.Linear(dim, num_classes))

    def forward(self, x):
        x = x.permute(0, 2, 1)
        x = x.unsqueeze(1)
        x = self.to_patch_embedding(x)
        b, n, _ = x.shape
        x += self.pos_embedding[:, :n]
        x = self.dropout(x)
        x = self.transformer(x)
        x = self.to_latent(x)
        x = x.reshape(-1, x.shape[2])
        return self.mlp_head(x)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class EvalItem:
    """One audio file paired with metadata for evaluation."""

    wav_path: Path
    label_path: Path
    sample_rate: float
    name: str
    noise: str = ""
    snr: int = 0

    @property
    def condition_key(self) -> str:
        return f"{self.noise}@{self.snr}dB"


# ---------------------------------------------------------------------------
# LibriVAD file discovery and label handling
# ---------------------------------------------------------------------------

def resolve_label_path(
    wav_path: Path,
    labels_root: Path,
    dataset: str,
    split: str,
) -> Path:
    """Map a LibriVAD noisy WAV path to its sample-level label .npy file."""
    stem = wav_path.stem
    parts = stem.split("_+_")[0].split("-")
    speaker = parts[0]
    chapter = parts[1]
    return labels_root / dataset / split / speaker / chapter / f"{stem}.npy"


def sample_labels_to_frame_labels(
    sample_labels: npt.NDArray[np.int16],
    frame_len_samples: int,
) -> npt.NDArray[np.int32]:
    """Downsample sample-level binary labels to frame-level via majority vote."""
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
) -> list[tuple[Path, str, int]]:
    """Find noisy WAV files in the LibriVAD Results tree."""
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
            for w in wavs:
                files.append((w, noise_name, snr_val))

    return files


def load_eval_items(
    librivad_root: Path,
    dataset: str,
    split: str,
    noises: list[str] | None,
    snrs: list[int] | None,
) -> list[EvalItem]:
    """Discover and validate LibriVAD files, returning lightweight EvalItems."""
    results_root = librivad_root / "Results"
    labels_root = librivad_root / "Files" / "Labels"

    discovered = discover_librivad_files(results_root, dataset, split, noises, snrs)
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

        info = sf.info(str(wav_path))
        items.append(EvalItem(
            wav_path=wav_path,
            label_path=label_path,
            sample_rate=info.samplerate,
            name=wav_path.stem,
            noise=noise_name,
            snr=snr_val,
        ))

    if skipped > 0:
        print(f"  Warning: skipped {skipped} files with missing labels.")

    return items


# ---------------------------------------------------------------------------
# Metric computation
# ---------------------------------------------------------------------------

def compute_auc(
    scores: npt.NDArray[np.float64],
    targets: npt.NDArray[np.int32],
) -> float | None:
    """Compute AUC-ROC, returning None if only one class is present."""
    if len(targets) == 0 or targets.min() == targets.max():
        return None
    return float(roc_auc_score(targets, scores))


def compute_fbeta(
    predictions: npt.NDArray[np.int32],
    targets: npt.NDArray[np.int32],
    beta: float = 2.0,
) -> dict[str, float]:
    """Compute F-beta score and supporting metrics."""
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


# ---------------------------------------------------------------------------
# miniDSP VAD evaluation
# ---------------------------------------------------------------------------

def eval_minidsp(
    items: list[EvalItem],
    frame_len_ms: float = MINIDSP_FRAME_MS,
) -> list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]]:
    """Run miniDSP VAD on all items, return (predictions, targets, scores) per item."""
    results = []
    for item in items:
        audio, sr = sf.read(str(item.wav_path), dtype="float64")
        if audio.ndim > 1:
            audio = audio.mean(axis=1)

        sample_labels = np.load(item.label_path).astype(np.int16)
        frame_len_samples = int(sr * frame_len_ms / 1000.0)

        min_samples = min(len(audio), len(sample_labels))
        audio = audio[:min_samples]
        sample_labels = sample_labels[:min_samples]

        targets = sample_labels_to_frame_labels(sample_labels, frame_len_samples)
        num_frames = len(targets)
        audio = audio[: num_frames * frame_len_samples]

        vad = VAD()
        decisions, scores, _ = vad.process(audio, sr, frame_len_samples)
        results.append((decisions, targets, scores))

    return results


# ---------------------------------------------------------------------------
# ViT-MFCC evaluation
# ---------------------------------------------------------------------------

def download_vit_model() -> Path:
    """Download the ViT-MFCC checkpoint from HuggingFace, return local path."""
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(
        repo_id=HF_REPO_ID,
        repo_type=HF_REPO_TYPE,
        filename=HF_MODEL_PATH,
        cache_dir=CACHE_DIR,
    )
    return Path(path)


def select_device() -> torch.device:
    """Pick best available PyTorch device: MPS > CUDA > CPU."""
    if hasattr(torch.backends, "mps") and torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def load_vit_model(checkpoint_path: Path, device: torch.device):
    """Instantiate the KWT model and load pre-trained weights."""
    model = KWT(
        input_res=VIT_INPUT_RES,
        patch_res=VIT_PATCH_RES,
        num_classes=VIT_NUM_CLASSES,
        dim=VIT_DIM,
        depth=VIT_DEPTH,
        heads=VIT_HEADS,
        mlp_dim=VIT_MLP_DIM,
        channels=VIT_CHANNELS,
        dim_head=VIT_DIM_HEAD,
        dropout=VIT_DROPOUT,
        emb_dropout=VIT_EMB_DROPOUT,
        pre_norm=VIT_PRE_NORM,
    )
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    model.to(device)
    model.eval()
    return model


def vit_infer_file(
    wav_path: Path,
    label_path: Path,
    model: torch.nn.Module,
    device: torch.device,
) -> tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]:
    """Run ViT-MFCC inference on a single file, return (predictions, targets, scores).

    Uses the reference MFCC extraction and sequence splitting from the
    LibriVAD Eval-ViT-MFCC code.  *scores* contains the per-frame speech
    probability from softmax.
    """
    from mfcc import mfcc, cmvn  # type: ignore[import-not-found]
    from auc_Vit import split_data_into_sequences  # type: ignore[import-not-found]

    # mfcc() expects a comma-separated string: "wav_path,label_path,feat_prefix"
    mfcc_input = f"{wav_path},{label_path},_unused"
    features, labels = mfcc(
        mfcc_input, MFCC_WINLEN, MFCC_OVRLEN, MFCC_PRE_COEF,
        MFCC_NFILTER, MFCC_NFFT, MFCC_NO_COEF,
    )
    features = cmvn(features)

    num_frames = features.shape[0]
    labels = np.asarray(labels).ravel()[:num_frames]

    X_sequences, y_sequences, last_feat, last_labels = split_data_into_sequences(
        features, labels.reshape(-1, 1), VIT_SEQUENCE_LENGTH,
    )
    remainder = num_frames % VIT_SEQUENCE_LENGTH

    chunks = list(X_sequences)
    if remainder > 0:
        chunks.append(last_feat)

    all_preds = []
    all_scores = []
    with torch.no_grad():
        for x_seq in chunks:
            x_t = torch.from_numpy(x_seq).unsqueeze(0).float().to(device)
            probs = torch.softmax(model(x_t), dim=1)
            all_preds.append(torch.argmax(probs, dim=1).cpu().numpy())
            all_scores.append(probs[:, 1].cpu().numpy())

    predictions = np.concatenate(all_preds)[:num_frames].astype(np.int32)
    scores_arr = np.concatenate(all_scores)[:num_frames].astype(np.float64)
    targets = labels[:num_frames].astype(np.int32)
    return predictions, targets, scores_arr


def eval_vit(
    items: list[EvalItem],
    model: torch.nn.Module,
    device: torch.device,
) -> list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]]:
    """Run ViT-MFCC on all items, return (predictions, targets, scores) per item."""
    results = []
    for item in items:
        preds, targets, scores = vit_infer_file(item.wav_path, item.label_path, model, device)
        results.append((preds, targets, scores))
    return results


# ---------------------------------------------------------------------------
# Aggregation and reporting
# ---------------------------------------------------------------------------

def aggregate_metrics(
    results: list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]],
    beta: float,
) -> dict[str, float | None]:
    """Aggregate per-file results into F-beta score and AUC-ROC.

    Reports two AUC variants:
      - ``auc_pooled``: computed on concatenated frames (micro-average)
      - ``auc_macro``: mean of per-file AUC values (macro-average, matching
        the methodology used in the LibriVAD paper Table 7)
    """
    preds = np.concatenate([r[0] for r in results])
    targets = np.concatenate([r[1] for r in results])
    scores = np.concatenate([r[2] for r in results])
    metrics = compute_fbeta(preds, targets, beta=beta)
    metrics["auc_pooled"] = compute_auc(scores, targets)

    # Macro-averaged AUC: per-file AUC then mean (skipping single-class files)
    per_file_aucs = [compute_auc(r[2], r[1]) for r in results]
    valid = [a for a in per_file_aucs if a is not None]
    metrics["auc_macro"] = float(np.mean(valid)) if valid else None
    return metrics



def _fmt_auc(val: float | None) -> str:
    """Format an AUC value, returning 'N/A' when undefined."""
    return "N/A" if val is None else f"{val:.4f}"


def print_overall(
    minidsp_metrics: dict[str, float | None],
    vit_metrics: dict[str, float | None],
    n_files: int,
    beta: float,
    minidsp_elapsed: float = 0.0,
    vit_elapsed: float = 0.0,
) -> None:
    """Print overall comparison table."""
    label = f"F{beta:g}"
    print(f"\nOverall ({n_files} files):")
    print(f"  {'System':<20s} {label:>8s}  {'Precision':>9s}  {'Recall':>6s}  {'AUC(m)':>8s}  {'AUC(p)':>8s}  {'Time':>8s}")
    print(f"  {'─' * 20} {'─' * 8}  {'─' * 9}  {'─' * 6}  {'─' * 8}  {'─' * 8}  {'─' * 8}")
    for name, m, elapsed in [
        ("miniDSP VAD", minidsp_metrics, minidsp_elapsed),
        ("ViT-MFCC (small)", vit_metrics, vit_elapsed),
    ]:
        print(
            f"  {name:<20s} {m['f_beta']:8.4f}  {m['precision']:9.4f}  {m['recall']:6.4f}"
            f"  {_fmt_auc(m.get('auc_macro')):>8s}  {_fmt_auc(m.get('auc_pooled')):>8s}  {elapsed:7.1f}s"
        )
    print(f"  AUC(m)=macro-averaged per file; AUC(p)=pooled across all frames")


def _bucket_items(items: list[EvalItem]) -> tuple[dict[str, list[int]], dict[int, list[int]]]:
    """Index items by noise type and by SNR level, returning {label: [indices]}."""
    by_noise: dict[str, list[int]] = defaultdict(list)
    by_snr: dict[int, list[int]] = defaultdict(list)
    for i, item in enumerate(items):
        by_noise[item.noise].append(i)
        by_snr[item.snr].append(i)
    return by_noise, by_snr


def _print_bucket_table(
    title: str,
    label_header: str,
    label_width: int,
    buckets: dict[str, list[int]],
    minidsp_results: list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]],
    vit_results: list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]],
    beta: float,
    format_label: callable,
) -> None:
    """Print a breakdown table grouped by the given buckets."""
    f_label = f"F{beta:g}"
    print(f"\n{title}:")
    print(f"  {label_header:<{label_width}s} {'miniDSP ' + f_label:>12s}  {'ViT ' + f_label:>12s}  {'miniDSP AUC':>12s}  {'ViT AUC':>12s}")
    print(f"  {'─' * label_width} {'─' * 12}  {'─' * 12}  {'─' * 12}  {'─' * 12}")
    for key in sorted(buckets):
        indices = buckets[key]
        md_m = aggregate_metrics([minidsp_results[i] for i in indices], beta)
        vt_m = aggregate_metrics([vit_results[i] for i in indices], beta)
        print(f"  {format_label(key):<{label_width}s} {md_m['f_beta']:12.4f}  {vt_m['f_beta']:12.4f}  {_fmt_auc(md_m.get('auc_macro')):>12s}  {_fmt_auc(vt_m.get('auc_macro')):>12s}")
    print(f"  AUC = macro-averaged per file")


def print_breakdown(
    items: list[EvalItem],
    minidsp_results: list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]],
    vit_results: list[tuple[npt.NDArray[np.int32], npt.NDArray[np.int32], npt.NDArray[np.float64]]],
    beta: float,
) -> None:
    """Print per-noise and per-SNR breakdown tables."""
    by_noise, by_snr = _bucket_items(items)

    _print_bucket_table(
        "Per Noise Type", "Noise", 20, by_noise,
        minidsp_results, vit_results, beta, str,
    )
    _print_bucket_table(
        "Per SNR", "SNR (dB)", 10, by_snr,
        minidsp_results, vit_results, beta, str,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Compare miniDSP VAD vs. ViT-MFCC on LibriVAD",
    )
    p.add_argument(
        "--librivad-root", type=Path, required=True,
        help="Path to LibriVAD project root (contains Results/ and Files/)",
    )
    p.add_argument(
        "--dataset", default="LibriSpeechConcat",
        choices=["LibriSpeech", "LibriSpeechConcat"],
        help="Dataset variant (default: LibriSpeechConcat)",
    )
    p.add_argument(
        "--split", default="test-clean",
        help="Dataset split (default: test-clean)",
    )
    p.add_argument(
        "--noises", nargs="+", default=None,
        help="Noise types to include (default: all)",
    )
    p.add_argument(
        "--snrs", nargs="+", type=int, default=None,
        help="SNR levels to include (default: all)",
    )
    p.add_argument(
        "--beta", type=float, default=2.0,
        help="F-beta parameter; >1 favors recall (default: 2.0)",
    )
    p.add_argument(
        "--breakdown", action="store_true",
        help="Print per-condition (noise/SNR) breakdown",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    librivad_root = args.librivad_root.expanduser().resolve()

    # Add LibriVAD Eval-ViT-MFCC to sys.path for model imports
    vit_code_dir = librivad_root / "Eval-ViT-MFCC"
    if not vit_code_dir.is_dir():
        print(f"Error: ViT-MFCC code not found at {vit_code_dir}", file=sys.stderr)
        print("Ensure the LibriVAD repo is at --librivad-root.", file=sys.stderr)
        sys.exit(1)
    sys.path.insert(0, str(vit_code_dir))

    print(f"=== VAD Comparison: {args.split} ({args.dataset}) ===")
    print(f"    beta={args.beta}")

    # --- Discover files ---
    print("\nDiscovering files...")
    items = load_eval_items(
        librivad_root, args.dataset, args.split, args.noises, args.snrs,
    )
    print(f"  Found {len(items)} files")

    # --- Download ViT model ---
    print("\nDownloading ViT-MFCC checkpoint (if needed)...")
    checkpoint_path = download_vit_model()
    device = select_device()
    print(f"  Checkpoint: {checkpoint_path}")
    print(f"  Device: {device}")

    # --- Load ViT model ---
    print("Loading ViT-MFCC model...")
    model = load_vit_model(checkpoint_path, device)

    # --- Evaluate miniDSP ---
    print(f"\nEvaluating miniDSP VAD ({MINIDSP_FRAME_MS}ms frames)...")
    t0 = time.perf_counter()
    minidsp_results = eval_minidsp(items)
    minidsp_elapsed = time.perf_counter() - t0
    minidsp_overall = aggregate_metrics(minidsp_results, args.beta)
    print(f"  F{args.beta:g}={minidsp_overall['f_beta']:.4f}  "
          f"P={minidsp_overall['precision']:.4f}  "
          f"R={minidsp_overall['recall']:.4f}  "
          f"AUC(m)={_fmt_auc(minidsp_overall.get('auc_macro'))}  "
          f"AUC(p)={_fmt_auc(minidsp_overall.get('auc_pooled'))}  "
          f"({minidsp_elapsed:.1f}s)")

    # --- Evaluate ViT-MFCC ---
    print(f"\nEvaluating ViT-MFCC ({MFCC_OVRLEN * 1000:.0f}ms frames)...")
    t0 = time.perf_counter()
    vit_results = eval_vit(items, model, device)
    vit_elapsed = time.perf_counter() - t0
    vit_overall = aggregate_metrics(vit_results, args.beta)
    print(f"  F{args.beta:g}={vit_overall['f_beta']:.4f}  "
          f"P={vit_overall['precision']:.4f}  "
          f"R={vit_overall['recall']:.4f}  "
          f"AUC(m)={_fmt_auc(vit_overall.get('auc_macro'))}  "
          f"AUC(p)={_fmt_auc(vit_overall.get('auc_pooled'))}  "
          f"({vit_elapsed:.1f}s)")

    # --- Results ---
    print_overall(minidsp_overall, vit_overall, len(items), args.beta,
                  minidsp_elapsed, vit_elapsed)

    if args.breakdown:
        print_breakdown(items, minidsp_results, vit_results, args.beta)


if __name__ == "__main__":
    main()
