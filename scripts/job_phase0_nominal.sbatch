#!/bin/bash
#SBATCH -J ceo2_phase0
#SBATCH -p eertekin
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH -o runs/phase0_nominal/%x_%j/slurm.out
#SBATCH -e runs/phase0_nominal/%x_%j/slurm.err
# #SBATCH -A <your_allocation>    # uncomment & set if SLURM asks for an account

set -Eeuo pipefail
IFS=$'\n\t'
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

# ---------- Resolve repo root & runner ----------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# If you ran sbatch from repo/scripts, repo root is its parent; otherwise it's SUBMIT_DIR itself.
if [[ -x "${SUBMIT_DIR}/scripts/run_nominal.sh" ]]; then
  REPO_DIR="${SUBMIT_DIR}"
elif [[ -x "$(dirname "${SUBMIT_DIR}")/scripts/run_nominal.sh" ]]; then
  REPO_DIR="$(dirname "${SUBMIT_DIR}")"
else
  echo "ERROR: Cannot locate scripts/run_nominal.sh starting from: ${SUBMIT_DIR}" >&2
  exit 127
fi

RUNNER="${REPO_DIR}/scripts/run_nominal.sh"
cd "${REPO_DIR}"
echo "SUBMIT_DIR: ${SUBMIT_DIR}"
echo "REPO_DIR:   ${REPO_DIR}"
echo "Runner:     ${RUNNER}"
ls -l "${RUNNER}" || { echo "ERROR: runner not found or not executable"; exit 127; }

# ---------- Toolchain (match your build) ----------
module purge
module load gcc/13.3.0
module load openmpi/5.0.1
module load fftw/3.3.10

# ---------- Python / Conda env (needed for PyYAML in the runner) ----------
# Use the site conda, same as on login shells.
source /sw/apps/anaconda3/2024.10/etc/profile.d/conda.sh
conda activate ceo2-shi

# Prove PyYAML is present (prints into logs)
python - <<'PY'
import sys, yaml
print(f"Python: {sys.version.split()[0]}")
print(f"PyYAML: {yaml.__version__}")
PY

# ---------- LAMMPS (shared) ----------
LAMMPS_ROOT=/u/gclopton/scratch/lammps/opt/lammps/stable_29Aug2024_update1
LMP_BIN="${LAMMPS_ROOT}/bin/lmp"
export LD_LIBRARY_PATH="${LAMMPS_ROOT}/lib64:${LAMMPS_ROOT}/lib:${LD_LIBRARY_PATH:-}"

# OpenMPI plugins & safe defaults (avoid PML hiccups)
# Keep threads low, let the module wire MPI correctly.
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
# (Optional) If you truly need to pin to TCP only, you can keep:
# export OMPI_MCA_btl=self,tcp
# Otherwise leave OMPI/OPAL env unset; the module provides the right paths.

# ---------- Dated run folder (always) ----------
STAMP="$(date +%b_%d_%H%M%S)"                         # e.g., Oct_22_164336
RUN_TAG="${STAMP}_${SLURM_JOB_NAME}_${SLURM_JOB_ID}"   # e.g., Oct_22_164336_ceo2_phase0_5480041
RUN_DIR="${REPO_DIR}/runs/phase0_nominal/${RUN_TAG}"
mkdir -p "${RUN_DIR}/logs" "${RUN_DIR}/dumps"

# Redirect ALL output into the dated folder
exec > >(tee -a "${RUN_DIR}/logs/job.out") 2> >(tee -a "${RUN_DIR}/logs/job.err" >&2)
echo "Run tag: ${RUN_TAG}"
echo "Nodes:   ${SLURM_NODELIST}"
echo "LAMMPS:  ${LMP_BIN}"

# ---------- Launch Phase 0 via runner ----------
# Force working dir so relative paths inside the runner work.
# Use bash -x for a crystal-clear trace in logs.
set +e
bash -x "${RUNNER}" --lmp "${LMP_BIN}" --run-dir "${RUN_DIR}"
RC=$?
set -e

# Copy SLURM’s own stdout/err into the dated folder for one-stop logs
SLURM_DIR="${REPO_DIR}/runs/phase0_nominal/${SLURM_JOB_NAME}_${SLURM_JOB_ID}"
for f in slurm.out slurm.err; do
  [[ -f "${SLURM_DIR}/${f}" ]] && cp -f "${SLURM_DIR}/${f}" "${RUN_DIR}/logs/${f}" || true
done

ln -sfn "${RUN_DIR}" "${REPO_DIR}/runs/phase0_nominal/LATEST"
echo "Done. See ${RUN_DIR}/logs and ${RUN_DIR}/dumps."
exit "${RC}"
