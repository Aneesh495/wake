#!/usr/bin/env python3
"""
Telemetry analysis and visualization script for WAKE CFD simulations.
Generates aerodynamic lift and drag histories, limit cycle phase portraits,
throughput curves, and FFT power spectral density (PSD) distributions.
"""

import sys
import os
import argparse
import csv

def parse_aero_csv(filepath):
    steps = []
    times = []
    cd = []
    cl = []
    cm = []

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            steps.append(int(row['step']))
            times.append(float(row['time']))
            cd.append(float(row['cd']))
            cl.append(float(row['cl']))
            cm.append(float(row['cm']))

    return {
        'step': steps,
        'time': times,
        'cd': cd,
        'cl': cl,
        'cm': cm
    }

def parse_perf_csv(filepath):
    steps = []
    mlups = []
    max_speed = []

    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            steps.append(int(row['step']))
            mlups.append(float(row['mlups']))
            max_speed.append(float(row['max_speed']))

    return {
        'step': steps,
        'mlups': mlups,
        'max_speed': max_speed
    }

def compute_fft(signal, dt=1.0):
    try:
        import numpy as np
        n = len(signal)
        if n < 8:
            return [], []
        # Remove DC component
        sig = np.array(signal) - np.mean(signal)
        fft_vals = np.fft.rfft(sig)
        freqs = np.fft.rfftfreq(n, d=dt)
        psd = np.abs(fft_vals) ** 2
        return freqs.tolist(), psd.tolist()
    except ImportError:
        return [], []

def plot_telemetry(aero_data, perf_data, output_path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("[WAKE PLOT] Matplotlib not installed; printing ASCII summary.")
        print(f"Total samples: {len(aero_data['step'])}")
        if aero_data['cd']:
            print(f"Mean Cd: {sum(aero_data['cd']) / len(aero_data['cd']):.4f}")
            print(f"Mean Cl: {sum(aero_data['cl']) / len(aero_data['cl']):.4f}")
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle("WAKE CFD Aerodynamic Telemetry and Analysis", fontsize=16, fontweight='bold')

    # Subplot 1: Force Coefficients vs Time
    ax1 = axes[0, 0]
    ax1.plot(aero_data['step'], aero_data['cd'], label=r'$C_D$ (Drag)', color='#d62728', lw=1.5)
    ax1.plot(aero_data['step'], aero_data['cl'], label=r'$C_L$ (Lift)', color='#1f77b4', lw=1.5)
    ax1.set_xlabel("Simulation Timestep")
    ax1.set_ylabel("Force Coefficient")
    ax1.set_title("Aerodynamic Force Histories")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')

    # Subplot 2: Phase Portrait (Cl vs Cd)
    ax2 = axes[0, 1]
    half_idx = len(aero_data['cd']) // 2
    ax2.plot(aero_data['cd'][half_idx:], aero_data['cl'][half_idx:], color='#9467bd', lw=1.5)
    ax2.set_xlabel(r'$C_D$ (Drag Coefficient)')
    ax2.set_ylabel(r'$C_L$ (Lift Coefficient)')
    ax2.set_title("Limit Cycle Phase Portrait")
    ax2.grid(True, alpha=0.3)

    # Subplot 3: Power Spectral Density
    ax3 = axes[1, 0]
    freqs, psd = compute_fft(aero_data['cl'][half_idx:])
    if freqs:
        ax3.plot(freqs, psd, color='#2ca02c', lw=1.5)
        ax3.set_xlabel("Frequency (1 / dt)")
        ax3.set_ylabel("Power Spectral Density")
        ax3.set_title("Lift Coefficient FFT Power Spectrum")
        ax3.set_yscale('log')
        ax3.grid(True, alpha=0.3)

    # Subplot 4: Computational Throughput (MLUPS)
    ax4 = axes[1, 1]
    if perf_data and perf_data['step']:
        ax4.plot(perf_data['step'], perf_data['mlups'], color='#ff7f0e', lw=1.5)
        ax4.set_xlabel("Simulation Timestep")
        ax4.set_ylabel("Throughput (MLUPS)")
        ax4.set_title("ISPC SPMD Execution Throughput")
        ax4.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    print(f"[WAKE PLOT] Saved telemetry plot to: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Plot WAKE CFD telemetry CSV data")
    parser.add_argument("--aero", type=str, required=True, help="Path to aerodynamics CSV")
    parser.add_argument("--perf", type=str, default=None, help="Path to performance CSV")
    parser.add_argument("--out", type=str, default="telemetry_plot.png", help="Output PNG path")
    args = parser.parse_args()

    if not os.path.exists(args.aero):
        print(f"Error: File not found: {args.aero}")
        sys.exit(1)

    aero = parse_aero_csv(args.aero)
    perf = parse_perf_csv(args.perf) if args.perf and os.path.exists(args.perf) else None
    plot_telemetry(aero, perf, args.out)

if __name__ == "__main__":
    main()
