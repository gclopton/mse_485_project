#!/usr/bin/env python3
import argparse, os, numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def fig_rho_radial(post, figs):
    rho = pd.read_csv(os.path.join(post, "rho_r_vs_frame.csv"))
    R   = pd.read_csv(os.path.join(post, "R_track_vs_frame.csv"))
    # columns: frame, then one column per r_center (Å) as strings
    frames = rho["frame"].to_numpy()
    rA = np.array([float(c) for c in rho.columns[1:]])
    rnm = rA/10.0
    # pick a few representative frames
    idxs = [0, max(1, len(frames)//3), max(2, 2*len(frames)//3), len(frames)-1]
    plt.figure()
    for f in idxs:
        y = rho.iloc[f,1:].to_numpy()
        plt.plot(rnm, y, label=f"frame {f}")
        Rt = R["R_track_A"].iloc[f]
        if not np.isnan(Rt):
            plt.axvline(Rt/10.0, ls="--", alpha=0.5)
    # draw threshold if available
    meta = os.path.join(post, "rho_meta.txt")
    if os.path.exists(meta):
        with open(meta) as fh:
            lines = {ln.split("=")[0].strip(): float(ln.split("=")[1]) for ln in fh if "=" in ln}
        thr = lines.get("threshold", None)
        if thr is not None:
            plt.axhline(thr, ls=":", alpha=0.7, label="0.9·rho0")
    plt.xlabel("r [nm]"); plt.ylabel("mass density [g/cm³]"); plt.legend()
    plt.tight_layout(); plt.savefig(os.path.join(figs,"rho_radial_overlays.png"), dpi=180)

def fig_ndef_time(post, figs):
    df = pd.read_csv(os.path.join(post, "ndef_vs_time.csv"))
    t  = df["t_ps"].to_numpy(); y = df.iloc[:,1].to_numpy()
    peak = y.max()
    # try to read t_half, t_1e
    t_half = t_1e = None
    meta = os.path.splitext(os.path.join(post,"ndef_vs_time.csv"))[0] + ".meta.txt"
    if os.path.exists(meta):
        with open(meta) as fh:
            for ln in fh:
                if ln.startswith("t_half_ps"): t_half = float(ln.split("=")[1])
                if ln.startswith("t_1e_ps"):   t_1e   = float(ln.split("=")[1])
    plt.figure()
    plt.plot(t, y, lw=2)
    plt.axhline(peak, ls=":", alpha=0.5)
    if t_half is not None: plt.axvline(t_half, ls="--", alpha=0.7, label=r"$t_{1/2}$")
    if t_1e  is not None:  plt.axvline(t_1e,  ls="--", alpha=0.5, label=r"$t_{1/e}$")
    plt.xlabel("t [ps]"); plt.ylabel(r"defect density [nm$^{-3}$]"); plt.legend()
    plt.tight_layout(); plt.savefig(os.path.join(figs,"defects_time_series.png"), dpi=180)

def fig_Tl_time(post, figs):
    df = pd.read_csv(os.path.join(post, "Tl_core.csv"))
    t  = df["t_ps"].to_numpy(); Tl = df["Tl_core_K"].to_numpy()
    plt.figure()
    plt.plot(t, Tl, lw=2)
    plt.xlabel("t [ps]"); plt.ylabel(r"$T_l$ core [K]")
    plt.tight_layout(); plt.savefig(os.path.join(figs,"Te_Tl_quench.png"), dpi=180)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--post", required=True, help="runs/.../post directory")
    ap.add_argument("--figs", required=True, help="runs/.../figures directory")
    ap.add_argument("--also_copy_to", default="", help="optional extra dir, e.g. analysis/phase1_minimal")
    args = ap.parse_args()
    os.makedirs(args.figs, exist_ok=True)
    fig_rho_radial(args.post, args.figs)
    fig_ndef_time(args.post, args.figs)
    fig_Tl_time(args.post, args.figs)
    if args.also_copy_to:
        os.makedirs(args.also_copy_to, exist_ok=True)
        for name in ["rho_radial_overlays.png","defects_time_series.png","Te_Tl_quench.png"]:
            src = os.path.join(args.figs, name)
            dst = os.path.join(args.also_copy_to, name)
            if os.path.exists(src):
                with open(src,'rb') as s, open(dst,'wb') as d: d.write(s.read())

if __name__ == "__main__":
    main()
