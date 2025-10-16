#!/usr/bin/env bash
# scripts/create_env.sh
set -Eeuo pipefail

SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
ENV_DIR="${ROOT}/env"
LOG_DIR="${ROOT}/logs"
mkdir -p "${LOG_DIR}"

MON="$(date -u +%b | tr '[:lower:]' '[:upper:]')" 
STAMP="${MON}_$(date -u +%d)_$(date -u +%H%M%SZ)"

LOG="${LOG_DIR}/create_env.${STAMP}.log"
exec > >(tee -a "${LOG}") 2>&1

echo "[env] repo root: ${ROOT}"
echo "[env] conda: $(command -v conda || true)"
conda --version || {
  echo "ERROR: conda not found"
  exit 1
}

YML_CANDIDATES=(
  "${ENV_DIR}/environment.yml"
  "${ENV_DIR}/environment.yaml"
)
YML=""
for f in "${YML_CANDIDATES[@]}"; do
  [[ -f "$f" ]] && {
    YML="$f"
    break
  }
done
[[ -n "${YML}" ]] || {
  echo "ERROR: no environment.yml or environment.yaml under ${ENV_DIR}"
  exit 2
}
echo "[env] environment file: ${YML}"


ENV_NAME="$(awk -F': *' '/^[[:space:]]*name:/ {print $2; exit}' "${YML}" || true)"
ENV_NAME="${ENV_NAME:-ceo2-shi}"
echo "[env] target env: ${ENV_NAME}"

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "[env] updating existing environment '${ENV_NAME}'"
  conda env update --name "${ENV_NAME}" --file "${YML}" --prune
else
  echo "[env] creating environment '${ENV_NAME}'"
  conda env create --file "${YML}"
fi

LOCK_FROM_HIST_TS="${ENV_DIR}/lock-from-history.${STAMP}.yml"
LOCK_EXPL_TS="${ENV_DIR}/explicit.${STAMP}.lock.txt"
LOCK_FROM_HIST="${ENV_DIR}/lock-from-history.yml"
LOCK_EXPL="${ENV_DIR}/explicit.lock.txt"

echo "[env] writing locks"
conda run -n "${ENV_NAME}" python -V
conda env export -n "${ENV_NAME}" --from-history >"${LOCK_FROM_HIST_TS}"
conda list -n "${ENV_NAME}" --explicit >"${LOCK_EXPL_TS}"

cp -f "${LOCK_FROM_HIST_TS}" "${LOCK_FROM_HIST}"
cp -f "${LOCK_EXPL_TS}" "${LOCK_EXPL}"

echo "[env] done. log -> ${LOG}"
echo "[env] locks -> ${LOCK_FROM_HIST_TS} , ${LOCK_EXPL_TS}"
