#!/usr/bin/env python3
"""
Plot the magnitude spectrum CSV produced by the magnitude_spectrum example.

Usage:
    python3 examples/plot_spectrum.py

Reads:  examples/magnitude_spectrum.csv
Writes: examples/magnitude_spectrum.png
"""

import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")  # non-interactive backend for PNG output
import matplotlib.pyplot as plt


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "magnitude_spectrum.csv")
    png_path = os.path.join(script_dir, "magnitude_spectrum.png")

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found. Run the magnitude_spectrum "
              f"example first.", file=sys.stderr)
        sys.exit(1)

    # Read the CSV data
    freqs = []
    mags = []
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            freqs.append(float(row["frequency_hz"]))
            mags.append(float(row["magnitude"]))

    # Create the figure with two subplots: linear and dB scale
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.suptitle("Magnitude Spectrum  (440 Hz + 1000 Hz + 2500 Hz)",
                 fontsize=14, fontweight="bold")

    # --- Top: linear magnitude ---
    ax1.plot(freqs, mags, color="#2563eb", linewidth=0.8)
    ax1.set_ylabel("Amplitude")
    ax1.set_title("Linear scale")
    ax1.grid(True, alpha=0.3)

    # Annotate the three peaks
    for target_freq in [440, 1000, 2500]:
        # Find the closest bin
        idx = min(range(len(freqs)), key=lambda i: abs(freqs[i] - target_freq))
        ax1.annotate(
            f"{target_freq} Hz",
            xy=(freqs[idx], mags[idx]),
            xytext=(freqs[idx] + 200, mags[idx] + 0.05),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="#666"),
            color="#333",
        )

    # --- Bottom: dB magnitude ---
    import numpy as np
    mags_arr = np.array(mags)
    # Floor at -120 dB to avoid log(0)
    mags_db = 20 * np.log10(np.maximum(mags_arr, 1e-6))
    ax2.plot(freqs, mags_db, color="#dc2626", linewidth=0.8)
    ax2.set_ylabel("Magnitude (dB)")
    ax2.set_xlabel("Frequency (Hz)")
    ax2.set_title("Logarithmic (dB) scale")
    ax2.set_ylim(bottom=-80)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    fig.savefig(png_path, dpi=150)
    print(f"Saved plot to {png_path}")


if __name__ == "__main__":
    main()
