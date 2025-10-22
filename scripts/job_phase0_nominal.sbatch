#!/bin/bash
#SBATCH --job-name=ceo2_phase0
#SBATCH --partition=mrsec           # adjust if your account uses a different queue
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1         # increase if you want MPI ranks > 1
#SBATCH --time=00:30:00
#SBATCH -o logs/slurm-%j.out
#SBATCH -e logs/slurm-%j.err

set -Eeuo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"

# --- Modules: match your cluster environment ---
module purge
module load gcc/13.3.0 || true
module load openmpi/5.0.1-gcc-13.3.0 || module load openmpi/5.0.1 || true
module load fftw3 || true

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

# --- LAMMPS binary: must include 'fix ttm/yaml' ---
# Set LMP_BIN to your YAML-enabled LAMMPS, e.g. from your group install.
LMP_BIN=${LMP_BIN:-/path/to/lmp-with-ttm-yaml}

# Optional sanity check inside the job environment
if ! "$LMP_BIN" -help 2>/dev/null | grep -qi "ttm/yaml"; then
  echo "ERROR: LAMMPS binary '$LMP_BIN' does not list 'ttm/yaml'." >&2
  echo "       Use a YAML-enabled build or adjust LMP_BIN." >&2
  exit 2
fi

# MPI ranks = SLURM_NTASKS (defaults to 1)
NP=${SLURM_NTASKS:-1}

# Launch Phase 0 nominal shot
srun scripts/run_nominal.sh --np "$NP" --lmp "$LMP_BIN"
