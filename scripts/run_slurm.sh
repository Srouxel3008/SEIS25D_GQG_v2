#!/usr/bin/env bash
# Sanitized multi-node Slurm example for the HPC-scale Marmousi VSP case.
# Build first with: make -C src -f Makefile.hpc

# ---------------- Slurm settings: adapt to your cluster ----------------

#SBATCH --job-name=seis25d_marmousi
#SBATCH --nodes=6
#SBATCH --ntasks=6
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=0
#SBATCH --time=48:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
##SBATCH --partition=YOUR_PARTITION
##SBATCH --account=YOUR_ACCOUNT

set -euo pipefail

# ---------------------- Settings users may change ----------------------

# Replace these generic module names with those provided by your cluster.
module purge
module load intel
module load intel-mpi
module load mkl
# module load mumps              # uncomment for a USE_MUMPS=1 build

# The script finds the repository from its own location. No absolute project
# path is required.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJ_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

# Select another published inversion by changing the MainInput filename. This
# explicit assignment can later be replaced by a SLURM_ARRAY_TASK_ID mapping.
INPUT_DIR="${PROJ_ROOT}/examples/Marmousi_VSP"
OUTPUT_ROOT="${PROJ_ROOT}/runs/Marmousi_VSP"
EXE="${PROJ_ROOT}/bin/Seis2D_GQG_FWI"

MainInput="${INPUT_DIR}/MainInput_C33.inp"
InitInput="${INPUT_DIR}/InitInput.inp"
InitInputQ="${INPUT_DIR}/InitInputQ.inp"
TrueInput="${INPUT_DIR}/TrueInput.inp"
TrueInputQ="${INPUT_DIR}/TrueInputQ.inp"

# -------------------------- Validation/setup ---------------------------

if [[ ! -x "${EXE}" ]]; then
    echo "ERROR: executable not found or not executable: ${EXE}" >&2
    echo "Build it before submission with: make -C src -f Makefile.hpc" >&2
    exit 2
fi

for input_file in "${MainInput}" "${InitInput}"; do
    if [[ ! -f "${input_file}" ]]; then
        echo "ERROR: mandatory input file not found: ${input_file}" >&2
        exit 2
    fi
done

for input_file in "${InitInputQ}" "${TrueInput}" "${TrueInputQ}"; do
    if [[ -n "${input_file}" && ! -f "${input_file}" ]]; then
        echo "ERROR: selected optional input file not found: ${input_file}" >&2
        exit 2
    fi
done

INPUT_FILES=(
    "${MainInput}"
    "${InitInput}"
    "${InitInputQ}"
    "${TrueInput}"
    "${TrueInputQ}"
)

job_id="${SLURM_JOB_ID:-manual}"
output_dir="${OUTPUT_ROOT}/slurm_${job_id}"
mkdir -p "${output_dir}"

configuration_summary="${output_dir}/configuration_summary.txt"
{
    echo "MAIN_INPUT_FILE=$(basename -- "${MainInput}")"
    echo "MAIN_INPUT_SHA256=$(sha256sum "${MainInput}" | awk '{print $1}')"
    echo "===== BEGIN MAIN INPUT ====="
    cat -- "${MainInput}"
    echo
    echo "===== END MAIN INPUT ====="
} > "${configuration_summary}"

submit_dir="${SLURM_SUBMIT_DIR:-${PWD}}"
slurm_stdout="${submit_dir}/slurm-${job_id}.out"
slurm_stderr="${submit_dir}/slurm-${job_id}.err"

copy_scheduler_logs() {
    while true; do
        [[ -f "${slurm_stdout}" ]] && cp -f "${slurm_stdout}" "${output_dir}/" || true
        [[ -f "${slurm_stderr}" ]] && cp -f "${slurm_stderr}" "${output_dir}/" || true
        sleep 300
    done
}

copy_final_logs() {
    if [[ -n "${LOG_COPY_PID:-}" ]]; then
        kill "${LOG_COPY_PID}" 2>/dev/null || true
        wait "${LOG_COPY_PID}" 2>/dev/null || true
    fi
    [[ -f "${slurm_stdout}" ]] && cp -f "${slurm_stdout}" "${output_dir}/" || true
    [[ -f "${slurm_stderr}" ]] && cp -f "${slurm_stderr}" "${output_dir}/" || true
}

copy_scheduler_logs &
LOG_COPY_PID=$!
trap copy_final_logs EXIT TERM INT

# One MPI rank is assigned to each node. Every rank uses the CPUs allocated to
# it through OpenMP/PARDISO.

# ----------------------- Threads and memory limit ----------------------

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${OMP_NUM_THREADS}}"
export NTHREAD_MAX="${NTHREAD_MAX:-${OMP_NUM_THREADS}}"
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE
export OMP_STACKSIZE="${OMP_STACKSIZE:-1G}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

if [[ -z "${MEM_MAX:-}" ]]; then
    mem_kib="$(awk '/^MemAvailable:/{print $2; exit}' /proc/meminfo)"
    MEM_MAX="$((mem_kib / 1024 / 1024))"
    (( MEM_MAX > 0 )) || MEM_MAX=1
fi
export MEM_MAX

# -------------------------- Environment record -------------------------

{
    echo "START_TIME=$(date --iso-8601=seconds)"
    echo "SLURM_JOB_ID=${job_id}"
    echo "SLURM_JOB_NODELIST=${SLURM_JOB_NODELIST:-unknown}"
    echo "SLURM_NTASKS=${SLURM_NTASKS:-1}"
    echo "SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-unknown}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
    echo "NTHREAD_MAX=${NTHREAD_MAX}"
    echo "MEM_MAX_GIB=${MEM_MAX}"
    echo "EXECUTABLE=${EXE}"
    echo "INPUT_DIR=${INPUT_DIR}"
    command -v mpiifx 2>/dev/null || true
    command -v srun 2>/dev/null || true
    lscpu 2>/dev/null || true
    free -h 2>/dev/null || true
    if command -v git >/dev/null 2>&1; then
        echo "GIT_COMMIT=$(git -C "${PROJ_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
        echo "GIT_STATUS=$(git -C "${PROJ_ROOT}" status --porcelain 2>/dev/null | wc -l) changed entries"
    fi
} > "${output_dir}/environment.txt"

# ------------------------------- Run -----------------------------------

cd "${output_dir}"
SECONDS=0
echo "===== SEIS25D_GQG Slurm run started: $(date) ====="
echo "Results: ${output_dir}"

set +e
srun --ntasks="${SLURM_NTASKS:-1}" "${EXE}" "${INPUT_FILES[@]}"
run_status=$?
set -e

printf 'Runtime: %02d:%02d:%02d\n' \
    "$((SECONDS / 3600))" "$(((SECONDS % 3600) / 60))" "$((SECONDS % 60))"
echo "Exit status: ${run_status}"

exit "${run_status}"
