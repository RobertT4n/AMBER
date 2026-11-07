#!/bin/bash
#SBATCH --job-name=CRAMqc
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --output=logs/CRAMqc_%j.out

#Date: 2026-04-27
#Purpose: CRAMqc — Ultima UG CRAM quality control report. Custom FastQC equivalent.
#         Runs 0.CRAMqc.py on all *.cram files in the current directory in parallel.
#         Outputs per-sample HTML (FastQC-style quality report) + TSV (Google Sheets ready).
#Usage:   sbatch 0.CRAMqc.sh              (all *.cram files in cwd)
#         bash   0.CRAMqc.sh              (interactive, uses GNU parallel)
#         python3 0.CRAMqc.py sample.cram --sample 200000 --amber   (single file, quick)
#         python3 /path/to/AMBER/Processing/1.Raw2Bam/0.CRAMqc.py \
#             442431-C1-CAGTCAGTTGCAGAT-CAGTCAGTTGCAGAT.1MIL.cram --sample 1000000 --threads 4 --amber
#INPUT:   *.cram files in current directory (unaligned Ultima CRAMs, no reference needed)
#OUTPUT (per sample):
#   *.CRAMqc.html  — FastQC-style HTML report (open in browser)
#   *.CRAMqc.tsv   — Summary stats table (paste → Google Sheets → split by Tab)
#REQUIRES: python3 + matplotlib + numpy + samtools in PATH + GNU parallel

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SAMPLE_READS=500000   # reads to sample per file (lower for speed, higher for accuracy)
PARALLEL_JOBS=4       # samples to run simultaneously (each uses THREADS_PER_JOB samtools threads)
THREADS_PER_JOB=2     # samtools threads per job (PARALLEL_JOBS × THREADS_PER_JOB ≤ total CPUs)

mkdir -p logs

_run_one() {
    local cram="$1"
    echo "[$(date '+%H:%M:%S')] START ${cram}"
    python3 "${SCRIPT_DIR}/0.CRAMqc.py" "${cram}" \
        --sample "${SAMPLE_READS}" \
        --threads "${THREADS_PER_JOB}" \
        --amber
    echo "[$(date '+%H:%M:%S')] DONE  ${cram}"
}
export -f _run_one
export SCRIPT_DIR SAMPLE_READS THREADS_PER_JOB

shopt -s nullglob; crams=(*.cram)
[[ ${#crams[@]} -eq 0 ]] && { echo "No *.cram files found."; exit 1; }

echo "════════════════════════════════════════════"
echo " CRAMqc  $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Samples : ${#crams[@]}  |  Parallel: ${PARALLEL_JOBS}  |  Reads: ${SAMPLE_READS}"
echo "════════════════════════════════════════════"

printf '%s\n' "${crams[@]}" | \
    parallel --no-notice -j "${PARALLEL_JOBS}" --joblog logs/CRAMqc_parallel.log _run_one {}

echo ""
echo "════════════════════════════════════════════"
echo " CRAMqc complete  $(date '+%H:%M:%S')"
echo "  HTML reports : *.CRAMqc.html"
echo "  Job log      : logs/CRAMqc_parallel.log"
echo "  Tip: paste *.CRAMqc.tsv into Google Sheets (split by Tab) to compare samples"
echo "════════════════════════════════════════════"
