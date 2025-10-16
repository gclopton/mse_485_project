#!/usr/bin/env bash
# scripts/build_lammps.sh
# Build and install LAMMPS with MPI + FFTW + KSPACE + EXTRA-FIX (TTM)
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
LOG_DIR="${ROOT}/logs"
mkdir -p "${LOG_DIR}"
TP_DIR="${ROOT}/.thirdparty"
mkdir -p "${TP_DIR}"

MON="$(date -u +%b | tr '[:lower:]' '[:upper:]')"
STAMP="${MON}_$(date -u +%d)_$(date -u +%H%M%SZ)"

LMP_REF="${LMP_REF:-stable}" 
PREFIX="${PREFIX:-$HOME/opt/lammps-ttm}" 
SRC_DIR="${TP_DIR}/lammps"
BUILD_DIR="${SRC_DIR}/build"

LOG="${LOG_DIR}/build_lammps.${STAMP}.log"
exec > >(tee -a "${LOG}") 2>&1

echo "[lammps] repo root: ${ROOT}"
echo "[lammps] preparing toolchain modules"
module purge || true
module load gcc || true
module load openmpi || true
module load cmake || true
module load fftw || module load fftw3 || true

echo "[lammps] compilers:"
which mpicc || true
which mpicxx || true
mpicc -v || true
mpicxx -v || true

echo "[lammps] fetching source into ${SRC_DIR}"
if [[ ! -d "${SRC_DIR}/.git" ]]; then
  git clone https://github.com/lammps/lammps.git "${SRC_DIR}"
fi
git -C "${SRC_DIR}" fetch --all --tags
git -C "${SRC_DIR}" checkout "${LMP_REF}"
echo "[lammps] source at: $(git -C "${SRC_DIR}" rev-parse --short HEAD)"

echo "[lammps] configuring CMake"
rm -rf "${BUILD_DIR}"
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

cmake ../cmake \
  -D CMAKE_INSTALL_PREFIX="${PREFIX}" \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=on \
  -D CMAKE_C_COMPILER=mpicc \
  -D CMAKE_CXX_COMPILER=mpicxx \
  -D FFT=FFTW3 \
  -D PKG_KSPACE=on \
  -D PKG_EXTRA-FIX=on

echo "[lammps] building"
cmake --build . -j

echo "[lammps] installing -> ${PREFIX}"
cmake --install .

export PATH="${PREFIX}/bin:${PATH}"
lmp -h | head -n 40

echo "[lammps] checking packages:"
lmp -h | grep -F "Installed packages" || true
lmp -h | grep -F "KSPACE" && echo "KSPACE present"
lmp -h | grep -F "EXTRA-FIX" && echo "EXTRA-FIX present"


echo "[lammps] validating parse of PPPM/ZBL/TTM"
bash "${ROOT}/scripts/validate_install.sh"

echo "[lammps] done -> ${LOG}"

