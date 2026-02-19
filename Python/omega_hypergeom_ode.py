#!/usr/bin/env python3
# ============================================================
# Finite-N derivatives of s_N(eps) from:
#   Omega_N(eps) = (2 eps)^((N-1)/2) / (1+2 eps) * 2F1(1, N-3/2; N/2; u)
#   u = 1/(1+2 eps)
#
# s_N(eps) = (1/N) log Omega_N(eps)
# s'_N(eps) and s''_N(eps) by centered finite difference
#
# Exports (NO HEADERS, purely numeric):
#   out/s1_raw_N{N}.dat
#   out/s1_mapped_N{N}.dat
#   out/s1_asymptotic.dat
#   out/s2_raw_N{N}.dat
#   out/s2_mapped_N{N}.dat
#   out/s2_asymptotic.dat
#   out/plot_s1_raw.(pdf/png)
#   out/plot_s1_mapped.(pdf/png)
#   out/plot_s2_raw.(pdf/png)
#   out/plot_s2_mapped.(pdf/png)
# ============================================================

import os
import math
from functools import lru_cache

import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp


# ---------------------------
# 1) Asymptotic curves (thermodynamic limit)
# ---------------------------
J = 1.0
ec = J / 2.0

def beta_asym(e: float) -> float:
    """
    Asymptotic beta(e) = s'(e) (piecewise for Berlin-Kac spherical model J=1).
    Consistent with your dbetaDe definition for s''.
    """
    if (-J/2.0) <= e <= ec:
        return 1.0 / (e + J/2.0)
    elif e >= ec:
        return 1.0 / (2.0 * e)
    else:
        return float("nan")

def dbetaDe(e: float) -> float:
    """Asymptotic s''(e) = d beta / d e (piecewise)."""
    if (-J/2.0) <= e <= ec:
        return -1.0 / (e + J/2.0)**2
    elif e >= ec:
        return -1.0 / (2.0 * e**2)
    else:
        return float("nan")


# ---------------------------
# 2) Precision + robust logOmega
# ---------------------------
def mp_dps_for_N(N: int) -> int:
    """
    Precision heuristic for LARGE N.
    You can increase these if you see NaNs/instabilities.
    """
    if N <= 4000:
        return 160
    if N <= 8000:
        return 240
    if N <= 16000:
        return 320
    if N <= 32000:
        return 420
    if N <= 64000:
        return 520
    if N <= 128000:
        return 650
    if N <= 256000:
        return 800
    if N <= 512000:
        return 950
    return 1100


def _eps_key(eps: float, ndigits: int = 14) -> float:
    """Cache key stabilization for floats."""
    return float(round(float(eps), ndigits))


def log_hyp2f1_robust(a, b, c, z, dps: int) -> mp.mpf:
    """
    Returns log( 2F1(a,b;c;z) ) robustly for large parameters.
    """
    mp.mp.dps = dps

    try:
        hyp = mp.hyp2f1(a, b, c, z, maxterms=200000)
        if hyp <= 0:
            hyp = abs(hyp)
        return mp.log(hyp)
    except mp.libmp.libhyper.NoConvergence:
        pass

    one_minus_z = mp.mpf(1) - z
    z2 = z / (z - mp.mpf(1))  # negative for z in (0,1)

    try:
        hyp2 = mp.hyp2f1(a, c - b, c, z2, maxterms=300000)
        if hyp2 <= 0:
            hyp2 = abs(hyp2)
        return (-a) * mp.log(one_minus_z) + mp.log(hyp2)
    except mp.libmp.libhyper.NoConvergence:
        hyp2 = mp.hyp2f1(a, c - b, c, z2, maxterms=800000)
        if hyp2 <= 0:
            hyp2 = abs(hyp2)
        return (-a) * mp.log(one_minus_z) + mp.log(hyp2)


@lru_cache(maxsize=None)
def logOmegaN_mp_cached(eps_key: float, N: int, dps: int) -> mp.mpf:
    """
    Cached log(Omega_N(eps)) with explicit dps in key.
    """
    if eps_key <= 0.0:
        return mp.ninf

    mp.mp.dps = dps

    eps_mp = mp.mpf(eps_key)
    u = mp.mpf(1) / (mp.mpf(1) + mp.mpf(2) * eps_mp)  # in (0,1) for eps>0

    a = mp.mpf(1)
    b = mp.mpf(N) - mp.mpf("1.5")   # N - 3/2
    c = mp.mpf(N) / mp.mpf(2)       # N/2

    loghyp = log_hyp2f1_robust(a, b, c, u, dps)

    log_pref = (mp.mpf(N - 1) / mp.mpf(2)) * mp.log(mp.mpf(2) * eps_mp) \
               - mp.log(mp.mpf(1) + mp.mpf(2) * eps_mp)

    return log_pref + loghyp


def sN(eps: float, N: int, dps: int) -> float:
    """s_N(eps) = (1/N) log Omega_N(eps)."""
    lo = logOmegaN_mp_cached(_eps_key(eps), N, dps)
    if lo == mp.ninf:
        return float("nan")
    return float(lo / mp.mpf(N))


def s1_centered_fd(eps: float, N: int, h: float, dps: int) -> float:
    """
    s'(eps) ~ [s(eps+h) - s(eps-h)] / (2 h)
    """
    sp = sN(eps + h, N, dps)
    sm = sN(eps - h, N, dps)
    if not (math.isfinite(sp) and math.isfinite(sm)):
        return float("nan")
    return (sp - sm) / (2.0 * h)


def s2_centered_fd(eps: float, N: int, h: float, dps: int) -> float:
    """
    s''(eps) ~ [s(eps+h) - 2 s(eps) + s(eps-h)] / h^2
    """
    sp = sN(eps + h, N, dps)
    s0 = sN(eps,     N, dps)
    sm = sN(eps - h, N, dps)
    if not (math.isfinite(sp) and math.isfinite(s0) and math.isfinite(sm)):
        return float("nan")
    return (sp - 2.0 * s0 + sm) / (h * h)


# ---------------------------
# 3) Mapping + I/O (no header)
# ---------------------------
def map_to_range(y: np.ndarray, yref_min: float, yref_max: float) -> np.ndarray:
    """Linear map of y into [yref_min, yref_max] using finite values only."""
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(y)
    if mask.sum() < 2:
        return y.copy()

    ymin = np.nanmin(y[mask])
    ymax = np.nanmax(y[mask])
    if not (np.isfinite(ymin) and np.isfinite(ymax)) or abs(ymax - ymin) < 1e-300:
        return y.copy()

    return (y - ymin) * (yref_max - yref_min) / (ymax - ymin) + yref_min


def save_dat_noheader(path: str, x: np.ndarray, y: np.ndarray):
    """Write 2 columns with NO header and NO comment lines."""
    arr = np.column_stack([x, y])
    np.savetxt(path, arr)  # pure numeric


# ---------------------------
# 4) Main
# ---------------------------
def main():
    # --- Ns ---
    Ns = [1000, 512000]

    # eps range
    eps_min = 0.18
    eps_max = 0.83
    npts = 260
    eps_grid = np.linspace(eps_min, eps_max, npts)

    # finite-difference steps (separati)
    h1 = 8.0e-4   # for s'
    h2 = 1.2e-3   # for s''  (your previous safe choice)

    out_dir = "out"
    os.makedirs(out_dir, exist_ok=True)

    # --- asymptotic curves ---
    y1_asym = np.array([beta_asym(float(e)) for e in eps_grid], dtype=float)
#    y2_asym = np.array([dbetaDe(float(e)) for e in eps_grid], dtype=float)

    save_dat_noheader(os.path.join(out_dir, "s1_asymptotic.dat"), eps_grid, y1_asym)
#    save_dat_noheader(os.path.join(out_dir, "s2_asymptotic.dat"), eps_grid, y2_asym)

    # reference ranges for mapping
    m1 = np.isfinite(y1_asym)
#    m2 = np.isfinite(y2_asym)
    y1ref_min, y1ref_max = float(np.nanmin(y1_asym[m1])), float(np.nanmax(y1_asym[m1]))
#    y2ref_min, y2ref_max = float(np.nanmin(y2_asym[m2])), float(np.nanmax(y2_asym[m2]))

    curves1_raw, curves1_map = {}, {}
    curves2_raw, curves2_map = {}, {}

    for N in Ns:
        dps = mp_dps_for_N(N)
        logOmegaN_mp_cached.cache_clear()

        print(f"\n=== N={N} | mp.dps={dps} | npts={npts} | h1={h1} | h2={h2} ===")

        y1_raw = np.empty_like(eps_grid, dtype=float)
#        y2_raw = np.empty_like(eps_grid, dtype=float)

        for i, eps in enumerate(eps_grid):
            e = float(eps)
            y1_raw[i] = s1_centered_fd(e, N, h1, dps)
            #y2_raw[i] = s2_centered_fd(e, N, h2, dps)

            if (i + 1) % max(1, npts // 10) == 0:
                print(f"  progress: {i+1}/{npts}")

        curves1_raw[N] = y1_raw
        #curves2_raw[N] = y2_raw

        y1_mapped = map_to_range(y1_raw, y1ref_min, y1ref_max)
        #y2_mapped = map_to_range(y2_raw, y2ref_min, y2ref_max)

        curves1_map[N] = y1_mapped
        #curves2_map[N] = y2_mapped

        save_dat_noheader(os.path.join(out_dir, f"s1_raw_N{N}.dat"), eps_grid, y1_raw)
        save_dat_noheader(os.path.join(out_dir, f"s1_mapped_N{N}.dat"), eps_grid, y1_mapped)

        #save_dat_noheader(os.path.join(out_dir, f"s2_raw_N{N}.dat"), eps_grid, y2_raw)
        #save_dat_noheader(os.path.join(out_dir, f"s2_mapped_N{N}.dat"), eps_grid, y2_mapped)

    # ---- plots s' ----
    plt.figure()
    plt.plot(eps_grid, y1_asym, linewidth=2.5, label="asymptotic")
    for N in Ns:
        plt.plot(eps_grid, curves1_raw[N], linewidth=1.6, label=f"N={N}")
    plt.xlabel("eps")
    plt.ylabel("s'_N(eps) (raw)")
    plt.title("Finite-N s'(eps) vs asymptotic (raw)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "plot_s1_raw.pdf"))
    plt.savefig(os.path.join(out_dir, "plot_s1_raw.png"), dpi=220)
    plt.close()

    plt.figure()
    plt.plot(eps_grid, y1_asym, linewidth=2.5, label="asymptotic")
    for N in Ns:
        plt.plot(eps_grid, curves1_map[N], linewidth=1.6, label=f"N={N} (mapped)")
    plt.xlabel("eps")
    plt.ylabel("s'_N(eps) (mapped)")
    plt.title("Finite-N s'(eps) vs asymptotic (mapped)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "plot_s1_mapped.pdf"))
    plt.savefig(os.path.join(out_dir, "plot_s1_mapped.png"), dpi=220)
    plt.close()

    # ---- plots s'' ----
    plt.figure()
    plt.plot(eps_grid, y2_asym, linewidth=2.5, label="asymptotic")
    for N in Ns:
        plt.plot(eps_grid, curves2_raw[N], linewidth=1.6, label=f"N={N}")
    plt.xlabel("eps")
    plt.ylabel("s''_N(eps) (raw)")
    plt.title("Finite-N s''(eps) vs asymptotic (raw)")
    plt.legend()
    plt.tight_layout()
#    plt.savefig(os.path.join(out_dir, "plot_s2_raw.pdf"))
#    plt.savefig(os.path.join(out_dir, "plot_s2_raw.png"), dpi=220)
    plt.close()

    # plt.figure()
    # plt.plot(eps_grid, y2_asym, linewidth=2.5, label="asymptotic")
    # for N in Ns:
    #     plt.plot(eps_grid, curves2_map[N], linewidth=1.6, label=f"N={N} (mapped)")
    # plt.xlabel("eps")
    # plt.ylabel("s''_N(eps) (mapped)")
    # plt.title("Finite-N s''(eps) vs asymptotic (mapped)")
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(os.path.join(out_dir, "plot_s2_mapped.pdf"))
    # plt.savefig(os.path.join(out_dir, "plot_s2_mapped.png"), dpi=220)
    # plt.close()

    print("\nDone.")
    print(f"Outputs written to ./{out_dir}/")
    print("  - s1_asymptotic.dat")
#    print("  - s2_asymptotic.dat")
    for N in Ns:
        print(f"  - s1_raw_N{N}.dat")
        print(f"  - s1_mapped_N{N}.dat")
#        print(f"  - s2_raw_N{N}.dat")
#        print(f"  - s2_mapped_N{N}.dat")
    print("  - plot_s1_raw.pdf / plot_s1_raw.png")
    print("  - plot_s1_mapped.pdf / plot_s1_mapped.png")
#    print("  - plot_s2_raw.pdf / plot_s2_raw.png")
#    print("  - plot_s2_mapped.pdf / plot_s2_mapped.png")


if __name__ == "__main__":
    main()
