# scripts/run_nominal.sh
#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
MON="$(date -u +%b | tr '[:lower:]' '[:upper:]')"
STAMP="${MON}_$(date -u +%d)_$(date -u +%H%M%SZ)"

TAG="${1:-phase0_nominal}"
RUNDIR="${ROOT}/runs/${TAG}/${STAMP}"
mkdir -p "${RUNDIR}"

module load fftw >/dev/null 2>&1 || true
export PATH="$HOME/opt/lammps-ttm/bin:$PATH"

VARFILE="${RUNDIR}/vars.lmp"
conda run -n ceo2-shi python "${ROOT}/scripts/yaml2lmpvars.py" "${ROOT}/config/physics.yaml" "${VARFILE}"

LMPOUT="${RUNDIR}/log.lammps"
SCREEN="${RUNDIR}/screen.txt"
INFILE="${ROOT}/inputs/lammps/in.main.lmp"

lmp -var VARFILE "${VARFILE}" -in "${INFILE}" -log "${LMPOUT}" -screen "${SCREEN}"
echo "Done → ${RUNDIR}"
