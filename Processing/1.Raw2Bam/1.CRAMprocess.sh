#!/bin/bash
#SBATCH --job-name=AMBER_CRAMprocess
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=1-1
#SBATCH --output=logs/CRAMprocess_%A_%a.out
#SBATCH --error=logs/CRAMprocess_%A_%a.err

#Date: 2026-04-27
#Purpose: AMBER pre-processing step 1 (Ultima UG). Unaligned CRAM → BCtagged BAM.
#         Filters reads lacking GCATGCGG bridge-ligation clamp, extracts SpatBC (15bp),
#         applies degenerate base filter, and appends SpatBC to QNAME for downstream UMI-tools.
#Usage:   sbatch /path/to/AMBER/Processing/1.Raw2Bam/1.CRAMprocess.sh   (from directory containing *.cram files)
#         bash   /path/to/AMBER/Processing/1.Raw2Bam/1.CRAMprocess.sh

#CHAIN:   Submit all CRAMs in parallel + auto-queue mapUG (run from your DATA directory):
#  [0] Load all modules on login node first:
#        module load bio/bwa-mem2/2.2.1-intel-compilers-2023.1.0 samtools cutadapt bedtools parallel
#  [1] Count CRAMs, submit as array (each task = 1 CRAM, all run simultaneously on cluster):
#        N=$(find . -maxdepth 1 -name "*.cram" | wc -l)
#        JID=$(sbatch --parsable --array=1-${N} /path/to/AMBER/Processing/1.Raw2Bam/1.CRAMprocess.sh)
#  [2] Immediately chain mapUG (also as array — waits for ALL CRAMprocess tasks to finish):
#        sbatch --array=1-${N} --dependency=afterok:$JID /path/to/AMBER/Processing/1.Raw2Bam/2.mapUG.sh [mm|hg|hgmm] bwamem2
#  Each CRAM processes independently with N_AWK_JOBS=8 awk workers internally.
#  afterok:$JID waits for ALL array tasks to succeed before starting mapUG array.

#         bash /path/to/AMBER/Processing/1.Raw2Bam/1.CRAMprocess.sh
#DEPENDS: samtools, GNU parallel  (load on login node: module load samtools parallel)
#INPUT:   Unaligned CRAM files (Ultima UG output, FLAG=4, CHROM=*, POS=0, MAPQ=0).
#         No reference genome needed — unaligned CRAMs contain no reference-compressed data.
#         Read structure: [14bp Tn5][~19bp Nextera ME ±1bp][15bp SpatBC][8bp clamp GCATGCGG][...]
#OUTPUT (per sample):
#   *.BCtagged.bam          — FINAL: unaligned BAM, SpatBC appended to QNAME as _<15bp>
#                             QNAME format: originalReadID_<SpatBC>, compatible with:
#                             UMI-tools --extract-umi-method=read_id --umi-separator='_'
#                             Custom tag hd:i:0 (exact clamp match) or hd:i:1 (1-mismatch).
#                             Split post-alignment: samtools view -d hd:0 / -d hd:1
#   *.CRAMprocess.log       — per-sample filter stats (clamp yield + SBC yield + overall)
#   logs/CRAMprocess_parallel.log  — job-level exit codes and timing
#NEXT:    Adapter-trim *.BCtagged.bam → align with mapUG.sh → 3.run_bamsplitUG.sh
#REQUIRES: samtools, GNU parallel

# ─── MODULES ─────────────────────────────────────────────────────────────────
module load samtools 2>/dev/null || true     # load on login node before sbatch; silenced in interactive

set -eo pipefail

# ─── CONFIG ───────────────────────────────────────────────────────────────────

# N_JOBS:     CRAMs processed simultaneously (outer parallelism).
#             Keep at 1 for large CRAMs (each CRAM uses all N_AWK_JOBS cores internally).
#             Increase only if CRAMs are small and you want them batched in parallel.
N_JOBS=1

# N_AWK_JOBS: awk workers per CRAM (inner parallelism via parallel --pipe).
#             parallel --pipe splits the SAM read stream into blocks and dispatches each
#             block exclusively to one worker — no read is ever processed by two workers.
#             Set to match --cpus-per-task (default: 8).
N_AWK_JOBS=8

# Bridge ligation clamp (8bp). Hard filter: reads without clamp are not valid AMBER molecules.
CLAMP="GCATGCGG"

# Clamp search window (1-based, inclusive).
# Tn5(14bp) + ME_min(18bp) + SpatBC(15bp) = 47bp minimum prefix → clamp cannot start <47.
# WIN_S=47: eliminates false early-exit on GCATGCGG appearing within the SpatBC region.
# WIN_E=65: accommodates +1bp ME jitter (peak at 49-50, spread 47-53, validated on 1M reads).
# CRITICAL: DO NOT hardcode clamp position — Nextera ME has insertion-length variation.
WIN_S=47; WIN_E=65

# Hamming distance tolerance for clamp matching (1 = absorbs single sequencing errors).
# Rationale: 8bp × 3 alternatives = 25 sequences match at hd≤1 out of 65,536 possible 8-mers.
# Exact match is ALWAYS preferred — 1-mismatch only accepted if no exact match in window.
MAX_MM=1

# SpatBC length (bp), extracted as the BC_LEN bases immediately before the clamp.
BC_LEN=15

# SpatBC degenerate design: BNNBNNDNNBNNBNN (B=not-A, D=not-C, N=any).
# B positions (1,4,10,13): reject A. D position (7): reject C.

# Input glob
INPUT_GLOB="*.cram"

# ──────────────────────────────────────────────────────────────────────────────

# Inner worker called by parallel --pipe for each SAM chunk.
# Receives reads-only SAM on stdin (no header lines — header handled separately).
# Writes one stats line (space-separated counts) to $1 for aggregation after all chunks finish.
_awk_chunk() {
    local statf="$1"
    awk -v clamp="${CLAMP}" -v clen=8 -v ws="${WIN_S}" -v we="${WIN_E}" \
        -v mm="${MAX_MM}" -v bclen="${BC_LEN}" -v statf="${statf}" \
        'BEGIN{ FS="\t"; OFS="\t"; n=0; ncp=0; nfail=0; npass=0; hd0=0; hd1=0 }
         {
           n++; seq=$10
           best_pos=0; best_hd=mm+1
           for(i=ws; i<=we; i++){
               obs=substr(seq,i,clen); hd=0
               for(j=1;j<=clen;j++) if(substr(obs,j,1)!=substr(clamp,j,1)) hd++
               if(hd<best_hd){ best_hd=hd; best_pos=i }
               if(hd==0) break
           }
           if(best_hd>mm || best_pos<=bclen) next
           ncp++; spatbc=substr(seq, best_pos-bclen, bclen)
           if(substr(spatbc,1,1)=="A" || substr(spatbc,4,1)=="A" || \
              substr(spatbc,7,1)=="C" || substr(spatbc,10,1)=="A" || \
              substr(spatbc,13,1)=="A"){ nfail++; next }
           $1=$1"_"spatbc; print $0 "\thd:i:" best_hd
           npass++; if(best_hd==0) hd0++; else hd1++
         }
         END{ print n, ncp, nfail, npass, hd0, hd1 > statf }'
}
export -f _awk_chunk

process_sample() {
    local cram="$1"
    local prefix="${cram%.cram}"
    local out_bam="${prefix}.BCtagged.bam"
    local report="${prefix}.CRAMprocess.log"
    local stats_dir; stats_dir=$(mktemp -d "${prefix}.stats_XXXXXX")

    echo "[START]  $cram  ($(date '+%H:%M:%S'))"

    # Fan reads across N_AWK_JOBS awk workers via parallel --pipe.
    # How it works: parallel reads the SAM stream, buffers 200MB blocks, sends each block
    # exclusively to one worker. No read is seen by more than one worker (parallel guarantees this).
    # Header (@ lines) extracted separately via samtools view -H and prepended by cat before
    # BAM conversion — workers only see read lines so header logic is not needed in awk.
    # -k preserves read order across workers (required for valid BAM).
    samtools view "${cram}" | \
    parallel --pipe --block 200M -k -j "${N_AWK_JOBS}" \
        "_awk_chunk ${stats_dir}/s_{#}.txt" | \
    cat <(samtools view -H "${cram}") - | \
    samtools view -b -o "${out_bam}"

    # Aggregate stats from all worker stat files into the final report
    awk -v label="${prefix}" \
        'BEGIN{ n=0; ncp=0; nfail=0; npass=0; hd0=0; hd1=0 }
         { n+=$1; ncp+=$2; nfail+=$3; npass+=$4; hd0+=$5; hd1+=$6 }
         END{
           printf "Sample: %s\n", label
           printf "  Total input reads      : %d\n", n
           printf "  [1] Clamp filter       : %d passed / %d total (%.1f%% retained, %.1f%% lost)\n", \
               ncp, n, (n>0?ncp/n*100:0), (n>0?(n-ncp)/n*100:0)
           printf "  [2] SBC degen filter   : %d passed / %d clamp-pass (%.1f%% retained, %.1f%% lost)\n", \
               npass, ncp, (ncp>0?npass/ncp*100:0), (ncp>0?nfail/ncp*100:0)
           printf "  [=] Overall yield      : %d / %d reads (%.1f%%)\n", \
               npass, n, (n>0?npass/n*100:0)
           printf "  [=] HD breakdown       : hd=0: %d, hd=1: %d (%.1f%% recovered via 1-mismatch)\n", \
               hd0, hd1, (npass>0?hd1/npass*100:0)
         }' "${stats_dir}"/s_*.txt > "${report}"

    rm -rf "${stats_dir}"
    cat "${report}"
    local n_reads; n_reads=$(samtools view -c "${out_bam}")
    echo "[DONE]   ${prefix}: ${n_reads} reads → $(basename "${out_bam}")  ($(date '+%H:%M:%S'))"
    echo ""
}

export -f process_sample
export CLAMP WIN_S WIN_E MAX_MM BC_LEN N_AWK_JOBS

mkdir -p logs

echo "════════════════════════════════════════════════════"
echo " 1.CRAMprocess.sh  $(date '+%Y-%m-%d %H:%M:%S')"
printf "  Job    : %-12s  Node: %s\n" "${SLURM_JOB_ID:-local}" "$(hostname -s)"
printf "  Config : CLAMP=%s  WIN=[%d,%d]  MAX_MM=%d  BC_LEN=%d  JOBS=%d  AWK_JOBS=%d\n" \
    "${CLAMP}" "${WIN_S}" "${WIN_E}" "${MAX_MM}" "${BC_LEN}" "${N_JOBS}" "${N_AWK_JOBS}"
printf "  samtools: %s\n" "$(samtools --version | head -1)"
printf "  Input  : %d CRAM files\n" "$(find . -maxdepth 1 -name "${INPUT_GLOB}" | wc -l)"
echo "════════════════════════════════════════════════════"
echo ""

# Array mode: each SLURM task processes exactly one CRAM (selected by task ID).
# Submit with: sbatch --array=1-N (see #CHAIN above). All tasks run simultaneously.
mapfile -t _CRAMS < <(find . -maxdepth 1 -name "${INPUT_GLOB}" | sort)
_TASK="${SLURM_ARRAY_TASK_ID:-1}"
_CRAM="${_CRAMS[$_TASK - 1]}"
[[ -z "${_CRAM}" ]] && { echo "ERROR: No CRAM for array task ${_TASK} (found ${#_CRAMS[@]} CRAMs)"; exit 1; }
process_sample "${_CRAM}"

# Append this sample's report to the shared summary log
echo "=== CRAMprocess summary ===" >> logs/CRAMprocess_summary.log
cat "${_CRAM%.cram}.CRAMprocess.log" >> logs/CRAMprocess_summary.log 2>/dev/null || true

echo ""
echo "════════════════════════════════════════════════════"
echo " 1.CRAMprocess.sh — COMPLETE  $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Output  : *.BCtagged.bam  (SpatBC in QNAME, hd:i tag)"
echo "  Reports : *.CRAMprocess.log  |  logs/CRAMprocess_summary.log"
echo "  Next    : 2.mapUG.sh"
echo "  Job log : logs/CRAMprocess_parallel.log"
echo "════════════════════════════════════════════════════"
