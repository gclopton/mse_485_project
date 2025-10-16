#!/usr/bin/env python3
"""

Usage
-----
python scripts/make_supercell.py --target-size 12 12 20

Behavior
--------
- Takes target box lengths in nm (Lx, Ly, Lz).
- Chooses integer unit-cell replications by ceiling along each axis.

Outputs
-------
- LAMMPS structure files are written to inputs/structure/data.ceo2_<Lx>x<Ly>x<Lz>nm.lmp
- Writes JSON manifest next to the data with the same basename.

Note: Default partial charges (+2.4 on Ce, -1.2 on O).
"""

from __future__ import annotations
import argparse, json, sys
from pathlib import Path
import numpy as np
from ase.io import read, write
from ase.build import make_supercell

# Repo paths
ROOT = Path(__file__).resolve().parents[1]
STRUCT_DIR = ROOT / "inputs" / "structure"
DEFAULT_POSCAR = STRUCT_DIR / "POSCAR_CeO2_unit.vasp"

# Masses (g/mol) — pinned for provenance
MASS_CE = 140.116
MASS_O  = 15.999

def parse_args():
    p = argparse.ArgumentParser(description="Build CeO2 supercell LAMMPS data from a unit POSCAR using a target size (nm).")
    p.add_argument("--target-size", nargs=3, type=float, metavar=("LX","LY","LZ"), required=True,
                   help="Target box lengths in nm. Script snaps to integer reps via ceiling (no strain).")
    p.add_argument("--poscar", type=Path, default=DEFAULT_POSCAR,
                   help=f"Path to unit-cell POSCAR (default: {DEFAULT_POSCAR})")
    args = p.parse_args()
    return args

def fmt_dim(x: float) -> str:
    """Pretty dimension label for filenames: 12 -> '12', 12.5 -> '12p5', 12.25 -> '12p25'."""
    if abs(x - round(x)) < 1e-9:
        return str(int(round(x)))
    # avoid dots in filenames
    s = f"{x:.3f}".rstrip("0").rstrip(".")
    return s.replace(".", "p")

def ensure_ce_o(atoms):
    species = sorted(set(atoms.get_chemical_symbols()))
    if species != ["Ce", "O"]:
        sys.stderr.write(f"[warn] Detected species {species}; expected ['Ce','O']. Continuing.\n")

def reorder_ce_first(atoms):
    """Deterministic type mapping: Ce->type 1, O->type 2 (stable order within species)."""
    syms = atoms.get_chemical_symbols()
    idx_ce = [i for i, s in enumerate(syms) if s == "Ce"]
    idx_o  = [i for i, s in enumerate(syms) if s == "O"]
    order_idx = np.array(idx_ce + idx_o, dtype=int)
    return atoms[order_idx]

def main():
    args = parse_args()
    if not args.poscar.exists():
        sys.exit(f"POSCAR not found: {args.poscar}")

    # Load unit cell
    atoms_uc = read(args.poscar.as_posix())
    ensure_ce_o(atoms_uc)

    # Target (nm) -> integer reps by ceiling (nm -> Å)
    targ_nm = np.array(list(map(float, args.target_size)), dtype=float)
    a_len, b_len, c_len = atoms_uc.cell.lengths()  # Å
    targ_A = targ_nm * 10.0
    nx = int(np.ceil(targ_A[0] / a_len))
    ny = int(np.ceil(targ_A[1] / b_len))
    nz = int(np.ceil(targ_A[2] / c_len))
    if min(nx, ny, nz) <= 0:
        sys.exit(f"Bad target {targ_nm} nm vs unit cell {a_len:.4f},{b_len:.4f},{c_len:.4f} Å.")

    # Build supercell
    rep = np.diag([nx, ny, nz])
    atoms_sc = make_supercell(atoms_uc, rep)

    # Deterministic type mapping and explicit masses/charges
    atoms_sc = reorder_ce_first(atoms_sc)
    syms_sc = np.array(atoms_sc.get_chemical_symbols())
    masses = np.where(syms_sc == "Ce", MASS_CE, MASS_O).astype(float)
    atoms_sc.set_masses(masses)
    # Simple, pinned partial charges (keep script minimal)
    charges = np.where(syms_sc == "Ce", 2.4, -1.2).astype(float)
    atoms_sc.set_initial_charges(charges)
    atoms_sc.set_pbc((True, True, True))
    post_L_A = np.array(atoms_sc.cell.lengths())  # Å (no padding)
    atoms_sc.wrap()

    # Auto-build output path from requested target size
    lx_lab, ly_lab, lz_lab = (fmt_dim(targ_nm[0]), fmt_dim(targ_nm[1]), fmt_dim(targ_nm[2]))
    out_name = f"data.ceo2_{lx_lab}x{ly_lab}x{lz_lab}nm.lmp"
    out_path = STRUCT_DIR / out_name
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Write LAMMPS data (atom_style charge)
    write(out_path.as_posix(), atoms_sc, format="lammps-data", atom_style="charge", units="metal")

    # Manifest (per-output, same basename)
    manifest_path = out_path.with_suffix(".json")
    n_atoms = int(len(atoms_sc))
    ce_count = int(np.count_nonzero(syms_sc == "Ce"))
    o_count  = int(np.count_nonzero(syms_sc == "O"))
    manifest = {
        "source_poscar": (str(args.poscar.relative_to(ROOT))
                          if str(args.poscar).startswith(str(ROOT)) else str(args.poscar)),
        "requested_target_nm": list(map(float, targ_nm)),
        "replicate": {"nx": nx, "ny": ny, "nz": nz},
        "box_lengths_nm": {
            "after_replication": list(map(lambda x: float(x)/10.0, post_L_A)),
        },
        "n_atoms": n_atoms,
        "composition_counts": {"Ce": ce_count, "O": o_count},
        "type_mapping": {"Ce": 1, "O": 2},
        "masses_g_per_mol": {"Ce": MASS_CE, "O": MASS_O},
        "charge_scheme": "partial",
        "charges_by_species": {"Ce": 2.4, "O": -1.2},
        "atom_style": "charge",
        "lammps_units": "metal",
        "output_data": (str(out_path.relative_to(ROOT))
                        if str(out_path).startswith(str(ROOT)) else str(out_path)),
    }
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Display Summary
    nm = post_L_A / 10.0
    print(f"[ok] Wrote {out_path}  "
          f"(box ~ {nm[0]:.3f} × {nm[1]:.3f} × {nm[2]:.3f} nm; "
          f"N={n_atoms}; replicate={nx}×{ny}×{nz}; charges=partial)")
    print(f"[ok] Manifest: {manifest_path}")

if __name__ == "__main__":
    main()
