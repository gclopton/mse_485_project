#!/usr/bin/env bash
set -Eeuo pipefail

TAG="${TAG:-stable_29Aug2024_update1}"
PREFIX="${PREFIX:-$PWD/opt/lammps/${TAG}}"
SRCROOT="${SRCROOT:-$PWD/.thirdparty}"
JOBS="${JOBS:-4}"
WANT_GPU="${WANT_GPU:-0}"
LMP_OUT="$PREFIX/bin/lmp"

echo "[build] TAG=${TAG}"
echo "[build] PREFIX=${PREFIX}"
echo "[build] SRCROOT=${SRCROOT}"
echo "[build] WANT_GPU=${WANT_GPU}"

# Mirror job's module set (safe to repeat)
module purge || true
module load gcc/13.3.0 || module load gcc || true
module load openmpi/5.0.1-gcc-13.3.0 || module load openmpi/5.0.1 || module load openmpi || true
module load cmake || true
module load fftw3 || module load fftw || true

# Temp dirs on $HOME to avoid locked /tmp
mkdir -p "$HOME/tmp" "$HOME/.conda/pkgs" "$(dirname "$PREFIX")" "$SRCROOT"
export TMPDIR="${TMPDIR:-$HOME/tmp}"

# Fetch LAMMPS
cd "$SRCROOT"
if [[ ! -d lammps ]]; then
  echo "[build] Cloning LAMMPS…"
  git clone https://github.com/lammps/lammps.git
fi
cd lammps
git fetch --all --tags
git checkout "$TAG"
echo "[build] LAMMPS commit: $(git rev-parse --short HEAD)"

# Configure (use ../cmake superbuild)
rm -rf build && mkdir build && cd build

CMAKE_FLAGS=(
  -D CMAKE_BUILD_TYPE=Release
  -D CMAKE_INSTALL_PREFIX="${PREFIX}"
  -D BUILD_MPI=ON
  -D BUILD_SHARED_LIBS=ON          # allows 'make install-python'
  -D CMAKE_CXX_STANDARD=17
  -D CMAKE_CXX_STANDARD_REQUIRED=ON
  # Packages needed for your workflow
  -D PKG_EXTRA-FIX=ON              # TTM fixes (ttm, ttm/grid, ttm/mod)
  -D PKG_EXTRA-DUMP=ON             # 'dump ... yaml'
  -D PKG_KSPACE=ON                 # PPPM etc.
  -D PKG_MOLECULE=ON
)

if [[ "$WANT_GPU" == "1" ]]; then
  echo "[build] GPU build requested (KOKKOS/CUDA) — only if your CC nodes have CUDA"
  CMAKE_FLAGS+=(
    -D PKG_KOKKOS=ON
    -D Kokkos_ENABLE_CUDA=ON
    -D Kokkos_ENABLE_OPENMP=ON
  )
else
  CMAKE_FLAGS+=(
    -D PKG_OPENMP=ON
    -D CMAKE_C_COMPILER=mpicc
    -D CMAKE_CXX_COMPILER=mpicxx
  )
fi

cmake "${CMAKE_FLAGS[@]}" ../cmake
cmake --build . -j"$JOBS"
cmake --install .

# Optional: install Python wheel into ACTIVE conda env (if any)
if command -v python >/dev/null 2>&1; then
  make install-python || true
else
  echo "[info] Skipping 'make install-python' (no Python/conda env active)."
fi

# Verify features (TTM + YAML); note: no 'fix ttm/yaml' in upstream
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,vader,tcp

if [[ ! -x "$LMP_OUT" ]]; then
  echo "ERROR: LAMMPS binary not found at $LMP_OUT"; exit 3
fi

echo "[verify] Banner:"
"$LMP_OUT" -h | head -n 5

echo "[verify] TTM fixes present?"
"$LMP_OUT" -h | sed -n '/^Fix styles/,/^Compute styles/p' | grep -E '(^| )(ttm|ttm/grid|ttm/mod)([ /]|$)' \
  || { echo "ERROR: TTM fixes missing (need PKG_EXTRA-FIX=ON)"; exit 2; }

echo "[verify] YAML thermo present?"
"$LMP_OUT" -h | sed -n '/^Thermo styles/,/^Compute styles/p' | grep -iq 'yaml' \
  || { echo "ERROR: thermo_style yaml missing (use a post-2022-03 tag)"; exit 2; }

echo "[verify] YAML dump present?"
"$LMP_OUT" -h | sed -n '/^Dump styles/,/^Compute styles/p' | grep -iq 'yaml' \
  || { echo "ERROR: dump yaml missing (need PKG_EXTRA-DUMP=ON)"; exit 2; }

# Tiny runtime smoke to prove YAML really writes
tmpd="$(mktemp -d)"; trap 'rm -rf "$tmpd"' EXIT
cat > "$tmpd/in.yaml_smoke" <<'EOF'
units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 4 0 4 0 4
create_box 1 box
create_atoms 1 box
pair_style lj/cut 2.5
pair_coeff * * 1.0 1.0 2.5
timestep 0.005
thermo 10
thermo_style yaml
dump d all yaml 20 dump.yaml id type x y z
dump_modify d thermo yes time yes units yes
run 40
EOF

"$LMP_OUT" -echo screen -log none -in "$tmpd/in.yaml_smoke" > "$tmpd/screen.out" 2>&1 || {
  echo "ERROR: runtime YAML smoke failed"; tail -n 50 "$tmpd/screen.out"; exit 4; }
grep -m1 '^thermo:' "$tmpd/dump.yaml" >/dev/null 2>&1 || echo "[info] YAML dump lacks embedded thermo (add 'dump_modify d thermo yes')."

echo "[done] LAMMPS installed: $LMP_OUT"
