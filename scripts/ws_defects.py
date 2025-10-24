#!/usr/bin/env python3
import argparse, os, glob
import numpy as np, pandas as pd

A_PER_NM = 10.0

def parse_single_frame(path):
    with open(path, "r") as f:
        while True:
            line = f.readline()
            if not line: break
            if not line.startswith("ITEM: TIMESTEP"): continue
            _ = int(f.readline().strip())
            assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
            natoms = int(f.readline().strip())
            assert f.readline().startswith("ITEM: BOX BOUNDS")
            xlo,xhi = map(float, f.readline().split()[:2])
            ylo,yhi = map(float, f.readline().split()[:2])
            zlo,zhi = map(float, f.readline().split()[:2])
            atoms_hdr = f.readline().strip()
            cols = atoms_hdr.split()[2:]
            idx = {c:i for i,c in enumerate(cols)}
            ids  = np.empty(natoms, dtype=np.int64)
            pos  = np.empty((natoms,3), dtype=float)
            for i in range(natoms):
                parts = f.readline().split()
                ids[i] = int(parts[idx["id"]])
                pos[i,0]= float(parts[idx["x"]]); pos[i,1]= float(parts[idx["y"]]); pos[i,2]= float(parts[idx["z"]])
            return (xlo,xhi,ylo,yhi,zlo,zhi), ids, pos
    raise RuntimeError(f"No frame in {path}")

def parse_sequence(pattern):
    frames=[]
    for path in sorted(glob.glob(pattern)):
        with open(path,"r") as f:
            while True:
                line=f.readline()
                if not line: break
                if not line.startswith("ITEM: TIMESTEP"): continue
                _ = int(f.readline().strip())
                assert f.readline().startswith("ITEM: NUMBER OF ATOMS")
                natoms = int(f.readline().strip())
                assert f.readline().startswith("ITEM: BOX BOUNDS")
                xlo,xhi = map(float, f.readline().split()[:2])
                ylo,yhi = map(float, f.readline().split()[:2])
                zlo,zhi = map(float, f.readline().split()[:2])
                atoms_hdr = f.readline().strip()
                cols = atoms_hdr.split()[2:]
                idx = {c:i for i,c in enumerate(cols)}
                ids  = np.empty(natoms, dtype=np.int64)
                pos  = np.empty((natoms,3), dtype=float)
                for i in range(natoms):
                    parts = f.readline().split()
                    ids[i] = int(parts[idx["id"]])
                    pos[i,0]= float(parts[idx["x"]]); pos[i,1]= float(parts[idx["y"]]); pos[i,2]= float(parts[idx["z"]])
                frames.append(((xlo,xhi,ylo,yhi,zlo,zhi), ids, pos))
    return frames

def minimum_image(d, L):  # wrap into [-L/2,L/2]
    return d - L*np.round(d/L)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pattern", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--frame_dt_ps", type=float, default=0.05)
    ap.add_argument("--rcyl_nm", type=float, default=8.0)
    ap.add_argument("--disp_thr_A", type=float, default=1.2)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    # Reference
    (xlo,xhi,ylo,yhi,zlo,zhi), ids_ref, pos_ref = parse_single_frame(args.ref)
    Lx,Ly,Lz = xhi-xlo, yhi-ylo, zhi-zlo
    cx, cy   = xlo+0.5*Lx, ylo+0.5*Ly
    rcyl_A   = args.rcyl_nm*A_PER_NM
    vol_nm3  = np.pi*(args.rcyl_nm**2)*(Lz/A_PER_NM)

    # Sequence
    seq = parse_sequence(args.pattern)
    id2i = {int(i):k for k,i in enumerate(ids_ref)}
    times = np.arange(len(seq))*args.frame_dt_ps
    ndef  = np.zeros(len(seq), dtype=float)

    for f, (box, ids, pos) in enumerate(seq):
        idxs = np.array([id2i[int(i)] for i in ids], dtype=int)
        d = pos - pos_ref[idxs]
        d[:,0] = minimum_image(d[:,0], Lx)
        d[:,1] = minimum_image(d[:,1], Ly)
        d[:,2] = minimum_image(d[:,2], Lz)
        disp = np.linalg.norm(d, axis=1)
        r = np.sqrt((pos[:,0]-cx)**2 + (pos[:,1]-cy)**2)
        mask = (r <= rcyl_A) & (disp > args.disp_thr_A)
        ndef[f] = mask.sum()/vol_nm3

    peak = float(ndef.max())
    i_peak = int(np.argmax(ndef))
    def first_cross(level):
        tail = np.where(ndef[i_peak:] <= level)[0]
        return (times[i_peak+tail[0]] if tail.size>0 else np.nan)
    t_half = first_cross(0.5*peak)
    t_1e   = first_cross(peak/np.e)

    df = pd.DataFrame({"t_ps": times, "defect_density_nm^-3": ndef})
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    with open(os.path.splitext(args.out)[0]+".meta.txt","w") as fh:
        fh.write(f"peak_defect_density_nm^-3 = {peak}\n")
        fh.write(f"t_half_ps = {t_half}\n")
        fh.write(f"t_1e_ps = {t_1e}\n")
        fh.write(f"rcyl_nm = {args.rcyl_nm}\n")
        fh.write(f"disp_thr_A = {args.disp_thr_A}\n")

if __name__ == "__main__":
    main()
