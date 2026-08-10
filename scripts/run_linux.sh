#!/usr/bin/env bash
# Local workstation launcher for WSL and native Linux.
#
# Usage:
#   ./scripts/run_linux.sh
#
# Optional runtime overrides:
#   OMP_NUM_THREADS=8
#   MKL_NUM_THREADS=8
#   MEM_MAX=64                  # available memory in GiB

set -euo pipefail

SECONDS=0

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

find_repository_root() {
    local candidate="${SCRIPT_DIR}"

    while [[ "${candidate}" != "/" ]]; do
        if [[ -d "${candidate}/examples/Small_TTI" && \
              -f "${candidate}/scripts/run_linux.sh" ]]; then
            printf '%s\n' "${candidate}"
            return 0
        fi
        candidate="$(dirname -- "${candidate}")"
    done

    return 1
}

if ! PROJ_ROOT="$(find_repository_root)"; then
    echo "ERROR: could not locate the repository root from ${SCRIPT_DIR}." >&2
    exit 2
fi

if (( $# != 0 )); then
    echo "ERROR: this launcher does not accept positional arguments." >&2
    echo "Usage: $0" >&2
    exit 2
fi

# ---------------------- Input files users may change -------------------

INPUT_DIR="${PROJ_ROOT}/examples/Small_TTI"
OUTPUT_ROOT="${PROJ_ROOT}/runs/Small_TTI"
EXE="${PROJ_ROOT}/bin/Seis2D_GQG_FWI"

MainInput="${INPUT_DIR}/MainInput.inp"
InitInput="${INPUT_DIR}/InitInput.inp"
InitInputQ="${INPUT_DIR}/InitInputQ.inp"
TrueInput="${INPUT_DIR}/TrueInput.inp"
TrueInputQ="${INPUT_DIR}/TrueInputQ.inp"

MPI_RANKS=1

if [[ -z "${OMP_NUM_THREADS:-}" ]]; then
    if [[ -r /proc/cpuinfo ]]; then
        OMP_NUM_THREADS="$(awk '/^processor[[:space:]]*:/{count++} END {print count+0}' /proc/cpuinfo)"
    elif command -v getconf >/dev/null 2>&1; then
        OMP_NUM_THREADS="$(getconf _NPROCESSORS_ONLN)"
    else
        OMP_NUM_THREADS=1
    fi
fi

if [[ ! "${OMP_NUM_THREADS}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: OMP_NUM_THREADS must be a positive integer." >&2
    exit 2
fi

if [[ -z "${MEM_MAX:-}" ]]; then
    if [[ -r /proc/meminfo ]]; then
        mem_kib="$(awk '/^MemAvailable:/{print $2; found=1; exit} /^MemTotal:/{total=$2} END {if (!found) print total}' /proc/meminfo)"
        MEM_MAX="$((mem_kib / 1024 / 1024))"
        (( MEM_MAX > 0 )) || MEM_MAX=1
    else
        MEM_MAX=1
    fi
fi

if [[ ! "${MEM_MAX}" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MEM_MAX must be a positive integer number of GiB." >&2
    exit 2
fi

export OMP_NUM_THREADS
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-${OMP_NUM_THREADS}}"
export NTHREAD_MAX="${NTHREAD_MAX:-${OMP_NUM_THREADS}}"
export MEM_MAX
export OMP_DYNAMIC=FALSE
export MKL_DYNAMIC=FALSE
export OMP_STACKSIZE="${OMP_STACKSIZE:-1G}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"

if ! command -v mpirun >/dev/null 2>&1; then
    echo "ERROR: mpirun was not found." >&2
    echo "Initialize the Linux Intel oneAPI and Intel MPI environment first." >&2
    exit 2
fi

if [[ ! -x "${EXE}" ]]; then
    echo "ERROR: executable not found or not executable: ${EXE}" >&2
    echo "Build the program before running this launcher." >&2
    exit 2
fi

# Pass all five roles explicitly. Until the Fortran argument interface is
# redesigned, skipping a middle positional file would shift every later role.
INPUT_FILES=(
    "${MainInput}"
    "${InitInput}"
    "${InitInputQ}"
    "${TrueInput}"
    "${TrueInputQ}"
)

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

timestamp="$(date +%Y%m%d_%H%M%S)"
run_start_kernel_time="$(date '+%Y-%m-%d %H:%M:%S')"
output_dir="${OUTPUT_ROOT}/output_${timestamp}"
mkdir -p "${output_dir}"

run_log="${output_dir}/run.log"
error_log="${output_dir}/error.log"
environment_log="${output_dir}/environment.txt"
configuration_summary="${output_dir}/configuration_summary.txt"

{
    echo "MAIN_INPUT_FILE=$(basename -- "${MainInput}")"
    echo "MAIN_INPUT_SHA256=$(sha256sum "${MainInput}" | awk '{print $1}')"
    echo "===== BEGIN MAIN INPUT ====="
    cat -- "${MainInput}"
    echo
    echo "===== END MAIN INPUT ====="
} > "${configuration_summary}"

{
    echo "===== Seis2D_GQG_FWI error log ====="
    echo "Start time: $(date)"
    echo "Normal output is recorded in run.log."
} > "${error_log}"

{
    echo "START_TIME=$(date --iso-8601=seconds)"
    echo "HOST=$(hostname)"
    echo "WORKSTATION_MODE=WSL_OR_LINUX"
    echo "EXECUTABLE=${EXE}"
    echo "INPUT_DIR=${INPUT_DIR}"
    echo "MPI_RANKS=${MPI_RANKS}"
    echo "OMP_NUM_THREADS=${OMP_NUM_THREADS}"
    echo "MKL_NUM_THREADS=${MKL_NUM_THREADS}"
    echo "NTHREAD_MAX=${NTHREAD_MAX}"
    echo "MEM_MAX_GIB=${MEM_MAX}"
    echo "OMP_DYNAMIC=${OMP_DYNAMIC}"
    echo "MKL_DYNAMIC=${MKL_DYNAMIC}"
    echo "OMP_STACKSIZE=${OMP_STACKSIZE}"
    echo "OMP_PROC_BIND=${OMP_PROC_BIND}"
    echo "OMP_PLACES=${OMP_PLACES}"
    echo "MPIRUN=$(command -v mpirun)"
    command -v mpiifx >/dev/null 2>&1 && echo "MPIIFX=$(command -v mpiifx)"
    command -v ifx >/dev/null 2>&1 && ifx --version | head -n 1
    mpirun --version 2>&1 | head -n 2 || true

    if command -v git >/dev/null 2>&1 && git -C "${PROJ_ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
        echo "GIT_COMMIT=$(git -C "${PROJ_ROOT}" rev-parse HEAD 2>/dev/null || echo unknown)"
        if git -C "${PROJ_ROOT}" diff --quiet --ignore-submodules -- 2>/dev/null && \
           git -C "${PROJ_ROOT}" diff --cached --quiet --ignore-submodules -- 2>/dev/null; then
            echo "GIT_STATUS=CLEAN"
        else
            echo "GIT_STATUS=DIRTY"
        fi
    else
        echo "GIT_COMMIT=nogit"
        echo "GIT_STATUS=unknown"
    fi
} > "${environment_log}"

echo "Starting Small_TTI example"
echo "Results: ${output_dir}"

cd "${output_dir}"
set +e
{
    echo "===== Seis2D_GQG_FWI workstation run ====="
    echo "Start time: $(date)"
    echo "MPI ranks: ${MPI_RANKS}"
    echo "OpenMP threads per rank: ${OMP_NUM_THREADS}"
    echo "MKL threads per rank: ${MKL_NUM_THREADS}"
    echo "------------------------------------------"
    mpirun -np "${MPI_RANKS}" "${EXE}" "${INPUT_FILES[@]}"
} > >(tee "${run_log}") 2> >(tee -a "${error_log}" >&2)
run_status=$?
set -e

duration=${SECONDS}
{
    echo "------------------------------------------"
    echo "End time: $(date)"
    echo "Exit status: ${run_status}"
    printf 'Total runtime: %02d:%02d:%02d\n' \
        $((duration / 3600)) $(((duration % 3600) / 60)) $((duration % 60))
} | tee -a "${run_log}"

# Record kernel out-of-memory events associated with this run when WSL/Linux
# permits ordinary users to read the kernel message buffer. This also catches
# messages such as "Killed process ... Seis2D_GQG_FWI".
set +e
kernel_messages="$(dmesg --since "${run_start_kernel_time}" 2>&1)"
dmesg_status=$?
set -e
if (( dmesg_status == 0 )); then
    oom_events="$(printf '%s\n' "${kernel_messages}" | \
        grep -Ei 'out of memory|oom-kill|killed process' || true)"
    if [[ -n "${oom_events}" ]]; then
        {
            echo
            echo "===== Kernel OOM evidence during this run ====="
            printf '%s\n' "${oom_events}"
        } | tee -a "${error_log}" >&2
    else
        echo "Kernel OOM check: no OOM event detected during this run." >> "${error_log}"
    fi
else
    {
        echo
        echo "Kernel OOM check unavailable: dmesg could not be read."
        printf '%s\n' "${kernel_messages}"
    } >> "${error_log}"
fi

if (( run_status != 0 )); then
    echo "Run failed. See ${run_log} and ${error_log}" | tee -a "${error_log}" >&2
    exit "${run_status}"
fi

echo "Run completed successfully: ${output_dir}"
