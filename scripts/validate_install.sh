#!/usr/bin/env bash
# scripts/validate_install.sh — Phase-0 validator (PPPM + buck/coul/long + TTM)
set -Eeuo pipefail

# 0) Find LAMMPS
LMP="$(command -v lmp || command -v lmp_mpi || command -v lmp_kokkos || command -v lmp_serial || true)"
[[ -z "$LMP" ]] && {
  echo "ERROR: no LAMMPS executable in PATH"
  exit 1
}

# 1) Confirm required packages in the banner (built-in)
HELP="$("$LMP" -h || true)"
grep -Fq "KSPACE" <<<"$HELP" || {
  echo "ERROR: KSPACE not listed"
  exit 2
}
grep -Fq "EXTRA-FIX" <<<"$HELP" || {
  echo "ERROR: EXTRA-FIX not listed"
  exit 3
}

# 2) Minimal input that initializes PPPM + buck/coul/long + fix ttm
#    Use a large box and ±1 charges so PPPM is happy; small real-space cutoff.
TMPDIR="$(mktemp -d)"
trap 'rm -rf "$TMPDIR"' EXIT
TMPIN="$TMPDIR/val.in"
TMPLOG="$TMPDIR/val.log"

cat >"$TMPIN" <<'EOF'
clear
units metal
atom_style charge
boundary p p p

# Large box so cutoffs/neighbors are well-behaved
region box block 0 50 0 50 0 50
create_box 1 box
create_atoms 1 single 10 10 10
create_atoms 1 single 40 40 40
mass 1 140.116
set atom 1 charge  1.0
set atom 2 charge -1.0

neighbor 2.0 bin
kspace_style pppm 1.0e-4

# Real-space Coulomb cutoff = 6 Å keeps lists modest
pair_style buck/coul/long 6.0
pair_coeff 1 1 1000.0 0.5 100.0   # A, rho, C (dummy; just for init)

group all type 1

# Minimal TTM (parameters arbitrary; parser/init only)
fix t all ttm 12345 1.0 1.0 10.0 0.1 0.0 2.0 1 1 1 set 300.0

run 0
EOF

if "$LMP" -echo none -log "$TMPLOG" -screen none -in "$TMPIN" >/dev/null 2>&1; then
  echo "OK: EXTRA-FIX & KSPACE present; PPPM/buck-coul-long/TTM parsed and initialized."
else
  echo "ERROR: parse/init failed. Last 80 lines:"
  tail -n 80 "$TMPLOG" || true
  exit 4
fi
