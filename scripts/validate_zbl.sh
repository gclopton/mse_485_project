#!/usr/bin/env bash
# scripts/validate_zbl.sh — smoke test for ZBL overlay with PPPM + Buckingham
set -Eeuo pipefail

LMP="$(command -v lmp || command -v lmp_mpi || command -v lmp_kokkos || command -v lmp_serial || true)"
[[ -z "$LMP" ]] && {
  echo "ERROR: no LAMMPS executable in PATH"
  exit 1
}

TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
IN="$TMPDIR/zbl_val.in"
LOG="$TMPDIR/zbl_val.log"

cat >"$IN" <<'EOF'
clear
units metal
atom_style charge
boundary p p p

# Large box + +/- charges so PPPM initializes cleanly
region box block 0 50 0 50 0 50
create_box 1 box
create_atoms 1 single 10 10 10
create_atoms 1 single 40 40 40
mass 1 140.116
set atom 1 charge  1.0
set atom 2 charge -1.0

neighbor 2.0 bin
kspace_style pppm 1.0e-4

# IMPORTANT: zbl needs inner/outer cutoffs here (e.g., 1.0 2.0 Å)
pair_style hybrid/overlay buck/coul/long 6.0 zbl 1.0 2.0

# ZBL needs atomic numbers in pair_coeff
pair_coeff 1 1 zbl 58 58                # Ce–Ce example
pair_coeff 1 1 buck/coul/long 1000.0 0.5 100.0   # dummy A rho C just to init

group all type 1
# Minimal TTM (parse/init only; numbers are placeholders)
fix t all ttm 12345 1.0 1.0 10.0 0.1 0.0 2.0 1 1 1 set 300.0

run 0
EOF

if "$LMP" -echo none -log "$LOG" -screen none -in "$IN" >/dev/null 2>&1; then
  echo "OK: ZBL overlay parsed and initialized with PPPM + buck/coul/long + TTM."
else
  echo "ERROR: ZBL parse/init failed. Last 80 lines:"
  tail -n 80 "$LOG" || true
  exit 4
fi
