#!/bin/bash
#SBATCH --job-name=AMBER_mapUG
#SBATCH --cpus-per-task=28
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --array=1-9          # ← SET THIS to the number of *.BCtagged.bam files in your directory
#SBATCH --output=logs/mapUG_%A_%a.out
#SBATCH --error=logs/mapUG_%A_%a.err

#Date: 2026-05-01
#Purpose: AMBER alignment step 2 (Ultima UG). BCtagged unaligned BAM → aligned + filtered BAMs.
#         Full pipe: samtools fastq | cutadapt (5'+3' trim, Q15) | aligner | samtools sort → raw BAM.
#         FASTQ is NEVER written to disk — exists only transiently in the pipe.
#         Then produces 3 filtered BAM variants for downstream bamsplit comparison.
#         SLURM array: each task processes ONE sample in parallel → all N samples finish in ~9h.
#Usage:
#  cd /your/data/directory   (must contain *.BCtagged.bam files)
#  N=$(ls *.BCtagged.bam | wc -l)
#  sbatch --array=1-${N} /path/to/AMBER/Processing/1.Raw2Bam/2.mapUG.sh mm bwamem2
#CHAIN:   Submit CRAMprocess + auto-queue mapUG in one go:
#  [0] module load bio/bwa-mem2/2.2.1-intel-compilers-2023.1.0 samtools cutadapt bedtools
#  [1] JID=$(sbatch --parsable .../1.CRAMprocess.sh)
#  [2] N=$(ls *.BCtagged.bam 2>/dev/null | wc -l || echo 9)
#      sbatch --dependency=afterok:${JID} --array=1-${N} .../2.mapUG.sh mm bwamem2
#         bash .../2.mapUG.sh mm bwamem2
#DEPENDS: bwamem2: module load bio/bwa-mem2/2.2.1-intel-compilers-2023.1.0 samtools cutadapt bedtools
#         bt2:     module load bowtie2 samtools cutadapt bedtools
#         hg   = hg38 only
#         mm   = mm10 only
#         hgmm = mm10+hg38 barnyard concat
#INPUT:   *.BCtagged.bam   (from 1.CRAMprocess.sh — unaligned BAM, SpatBC in QNAME as _<15bp>)
#OUTPUT (per sample):
#   *.{aligner}.raw.bam          — aligned, coordinate-sorted
#   *.{aligner}.q20.fullRmsk.bam — MAPQ≥20 + full milliDiv<140 RepeatMasker filter [PRIMARY]
#   *.{aligner}.q0.fullRmsk.bam  — no MAPQ filter + full RepeatMasker filter
#   *.{aligner}.q0.noRmsk.bam    — no MAPQ + no repeat filter
#   *.{aligner}.flagstat.log / *.cutadapt.log / *.align.log / *.{aligner}.stats.log
#   logs/mapUG_row_{prefix}.tsv  — per-sample stats row (cat all → logs/mapUG_summary.tsv)
#NEXT:    Feed filtered BAMs into 3.run_bamsplitUG.sh
#
# ── Aligner notes ─────────────────────────────────────────────────────────────
# bwamem2 (default): handles Ultima 300bp natively, propagates hd:i tag via -T hd + -C.
# bt2 (comparison only): designed for ≤250bp Illumina; hd:i tag NOT propagated.
# ─────────────────────────────────────────────────────────────────────────────
#
# ── cutadapt trimming notes ───────────────────────────────────────────────────
# TWO-PASS strategy (required to cleanly remove the full 5' non-genomic prefix):
#   Pass 1: -g GTATAAGAGACAG removes Tn5+ME (using last 13nt; first 6nt AGATGT unreliable
#           due to Ultima homopolymer indels). Untrimmed reads pass through unchanged.
#   Pass 2: -g GCATGCGG removes SpatBC+clamp → clean insert only.
# 3' trim: -a CTGTCTCTTATACACAT (17nt Nextera N5 rev-comp; pos 17-18 unreliable in Ultima)
# Quality: -q 15 (Ultima Q15 ≈ Illumina Q20 for effective error rate)
# ─────────────────────────────────────────────────────────────────────────────

set -eo pipefail

# ─── ARGS ─────────────────────────────────────────────────────────────────────
GENOME="${1:?ERROR: Usage: sbatch --array=1-N 2.mapUG.sh [hg|mm|hgmm] [bwamem2|bt2]}"
ALIGNER="${2:-bwamem2}"
[[ "${ALIGNER}" == "bwamem2" || "${ALIGNER}" == "bt2" ]] || \
    { echo "ERROR: Unknown aligner '${ALIGNER}'. Use: bwamem2 | bt2"; exit 1; }

# ─── MODULES (load on login node BEFORE sbatch) ───────────────────────────────
module load samtools 2>/dev/null || true
module load cutadapt 2>/dev/null || true
module load bedtools 2>/dev/null || true
module load bio/bwa-mem2/2.2.1-intel-compilers-2023.1.0 2>/dev/null || true

# The bare bwa-mem2 dispatcher can pick an AVX512 build the node cannot run (job dies silently).
# Use the AVX2 binary when present; override with BWAMEM2_BIN=... if needed.
BWAMEM2_BIN="${BWAMEM2_BIN:-$(command -v bwa-mem2.avx2 || command -v bwa-mem2)}"

# ─── THREAD CONFIG ────────────────────────────────────────────────────────────
T_FASTQ=2; T_CUTADAPT=4; T_ALIGN=20; T_SORT=4

# ─── GENOME CONFIG ────────────────────────────────────────────────────────────
# Reference paths come from environment variables (export before sbatch). Build instructions: resources/README.md
#   HG38_FASTA / MM10_FASTA / HGMM_FASTA : genome FASTA (bwa-mem2 index files must sit next to it)
#   HG38_RMSK  / MM10_RMSK  / HGMM_RMSK  : RepeatMasker bed (milliDiv < 140), built by resources/build_repeatmasker_bed.sh
case "${GENOME}" in
  hg)
    REF_BWAMEM2="${HG38_FASTA:?Set HG38_FASTA (see resources/README.md)}"
    REF_BT2="${HG38_BT2:-}"   # bowtie2 index prefix (optional; only for --aligner bt2 comparison)
    RMSK="${HG38_RMSK:?Set HG38_RMSK (see resources/README.md)}"
    KEEP_CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
    ;;
  mm)
    REF_BWAMEM2="${MM10_FASTA:?Set MM10_FASTA (see resources/README.md)}"
    REF_BT2="${MM10_BT2:-}"
    RMSK="${MM10_RMSK:?Set MM10_RMSK (see resources/README.md)}"
    KEEP_CHRS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY"
    ;;
  hgmm)
    REF_BWAMEM2="${HGMM_FASTA:?Set HGMM_FASTA (see resources/README.md)}"
    REF_BT2="${HGMM_BT2:-}"
    RMSK="${HGMM_RMSK:?Set HGMM_RMSK (see resources/README.md)}"
    KEEP_CHRS="hg_chr1 hg_chr2 hg_chr3 hg_chr4 hg_chr5 hg_chr6 hg_chr7 hg_chr8 hg_chr9 hg_chr10 hg_chr11 hg_chr12 hg_chr13 hg_chr14 hg_chr15 hg_chr16 hg_chr17 hg_chr18 hg_chr19 hg_chr20 hg_chr21 hg_chr22 hg_chrX hg_chrY mm_chr1 mm_chr2 mm_chr3 mm_chr4 mm_chr5 mm_chr6 mm_chr7 mm_chr8 mm_chr9 mm_chr10 mm_chr11 mm_chr12 mm_chr13 mm_chr14 mm_chr15 mm_chr16 mm_chr17 mm_chr18 mm_chr19 mm_chrX mm_chrY"
    ;;
  *) echo "ERROR: Unknown genome '${GENOME}'. Use: hg | mm | hgmm"; exit 1 ;;
esac

# ─── PICK THIS TASK'S SAMPLE ─────────────────────────────────────────────────
mkdir -p logs
INPUT_FILES=(*.BCtagged.bam)
N=${#INPUT_FILES[@]}
TASK_ID=${SLURM_ARRAY_TASK_ID:-1}   # default to 1 for local dry-runs
if [[ ${TASK_ID} -gt ${N} ]]; then
    echo "Array task ${TASK_ID} > ${N} BCtagged.bam files — nothing to do."; exit 0
fi
BC_BAM="${INPUT_FILES[$TASK_ID-1]}"
PREFIX="${BC_BAM%.BCtagged.bam}"
RAW_BAM="${PREFIX}.${ALIGNER}.raw.bam"

# ─── HEADER ──────────────────────────────────────────────────────────────────
echo "════════════════════════════════════════════════════"
echo " 2.mapUG.sh  $(date '+%Y-%m-%d %H:%M:%S')"
printf "  Job     : %-12s  Array task: %s/%s\n" "${SLURM_JOB_ID:-local}" "${TASK_ID}" "${N}"
printf "  Node    : %s\n" "$(hostname -s)"
printf "  Genome  : %s   Aligner: %s\n" "${GENOME}" "${ALIGNER}"
[[ "${ALIGNER}" == "bwamem2" ]] && printf "  REF     : %s\n" "${REF_BWAMEM2}" || printf "  REF     : %s\n" "${REF_BT2}"
printf "  RMSK    : %s\n" "${RMSK}"
printf "  Threads : fastq=%d  cutadapt=%d  align=%d  sort=%d\n" "${T_FASTQ}" "${T_CUTADAPT}" "${T_ALIGN}" "${T_SORT}"
printf "  samtools: %s\n" "$(samtools --version | head -1)"
printf "  cutadapt: %s\n" "$(cutadapt --version 2>/dev/null | head -1 || echo 'unknown')"
[[ "${ALIGNER}" == "bwamem2" ]] && \
    printf "  bwa-mem2: %s\n" "$("${BWAMEM2_BIN}" version 2>/dev/null | head -1 || echo 'unknown')" || \
    printf "  bowtie2 : %s\n" "$(bowtie2 --version 2>/dev/null | head -1 || echo 'unknown')"
printf "  Input   : %s\n" "${BC_BAM}"
echo "════════════════════════════════════════════════════"
echo ""

# ─── SKIP if valid raw BAM already exists ─────────────────────────────────────
if [[ -f "${RAW_BAM}" ]] && samtools quickcheck "${RAW_BAM}" 2>/dev/null; then
    echo "[SKIP]   ${PREFIX} — valid ${RAW_BAM} exists"; exit 0
fi

# ─── ALIGN ───────────────────────────────────────────────────────────────────
echo "[START]  ${BC_BAM}  ($(date '+%H:%M:%S'))"
if [[ "${ALIGNER}" == "bwamem2" ]]; then
    samtools fastq -@ "${T_FASTQ}" -T hd "${BC_BAM}" | \
    cutadapt -g GTATAAGAGACAG -e 0.18 -j 2 - 2>/dev/null | \
    cutadapt -g GCATGCGG -a CTGTCTCTTATACACAT -e 0.18 -q 15 --minimum-length 20 \
        -j "${T_CUTADAPT}" - 2>"${PREFIX}.cutadapt.log" | \
    "${BWAMEM2_BIN}" mem -t "${T_ALIGN}" -C "${REF_BWAMEM2}" - 2>"${PREFIX}.align.log" | \
    samtools sort -@ "${T_SORT}" -o "${RAW_BAM}"
else
    samtools fastq -@ "${T_FASTQ}" "${BC_BAM}" | \
    cutadapt -g GTATAAGAGACAG -e 0.18 -j 2 - 2>/dev/null | \
    cutadapt -g GCATGCGG -a CTGTCTCTTATACACAT -e 0.18 -q 15 --minimum-length 20 \
        -j "${T_CUTADAPT}" - 2>"${PREFIX}.cutadapt.log" | \
    bowtie2 -x "${REF_BT2}" -U - -p "${T_ALIGN}" 2>"${PREFIX}.align.log" | \
    samtools sort -@ "${T_SORT}" -o "${RAW_BAM}"
fi
samtools index "${RAW_BAM}"

# ─── RAW BAM STATS ───────────────────────────────────────────────────────────
samtools flagstat "${RAW_BAM}" > "${PREFIX}.${ALIGNER}.flagstat.log"
_mapped=$(samtools view -c -F4 "${RAW_BAM}")
_map_pct=$(grep "mapped (" "${PREFIX}.${ALIGNER}.flagstat.log" | head -1 | grep -oE '[0-9]+\.[0-9]+%' || echo "?%")
_trim_out=$(grep "Reads written" "${PREFIX}.cutadapt.log" | sed 's/.*: *//' | tr -d ',' || echo "?")
echo "[TRIM]   ${PREFIX}: cutadapt passed ${_trim_out}"
echo "[ALIGN]  ${PREFIX}: ${_mapped} mapped reads (${_map_pct}) → $(basename ${RAW_BAM})"

samtools stats "${RAW_BAM}" | grep "^SN" | \
    grep -E "raw total|reads mapped|average length|bases mapped \(cigar\)|error rate|average quality" \
    > "${PREFIX}.${ALIGNER}.stats.log"

_hd0="n/a"; _hd1="n/a"; _pct_hd1="n/a"
if [[ "${ALIGNER}" == "bwamem2" ]]; then
    _hd0=$(samtools view -c -d hd:0 "${RAW_BAM}")
    _hd1=$(samtools view -c -d hd:1 "${RAW_BAM}")
    _pct_hd1=$(awk "BEGIN{t=${_hd0}+${_hd1}; printf \"%.1f\", (t>0?${_hd1}/t*100:0)}")
    echo "[HD-raw] ${PREFIX}: hd=0: ${_hd0}  hd=1: ${_hd1}  (${_pct_hd1}% via 1-mismatch)"
fi

# ─── FILTERS A/B/C (run in parallel — all read from RAW_BAM simultaneously) ───
BAM_A="${PREFIX}.${ALIGNER}.q20.fullRmsk.bam"
BAM_B="${PREFIX}.${ALIGNER}.q0.fullRmsk.bam"
BAM_C="${PREFIX}.${ALIGNER}.q0.noRmsk.bam"

samtools view -b -F 2308 -q 20 -@ 4 "${RAW_BAM}" ${KEEP_CHRS} | \
    bedtools intersect -v -abam stdin -b "${RMSK}" > "${BAM_A}" &  PID_A=$!
samtools view -b -F 2308 -@ 4 "${RAW_BAM}" ${KEEP_CHRS} | \
    bedtools intersect -v -abam stdin -b "${RMSK}" > "${BAM_B}" &  PID_B=$!
samtools view -b -F 2308 -@ 4 "${RAW_BAM}" ${KEEP_CHRS} > "${BAM_C}" &  PID_C=$!

wait $PID_A || { echo "ERROR: FILT-A failed"; exit 1; }
wait $PID_B || { echo "ERROR: FILT-B failed"; exit 1; }
wait $PID_C || { echo "ERROR: FILT-C failed"; exit 1; }

samtools index "${BAM_A}" & samtools index "${BAM_B}" & samtools index "${BAM_C}" & wait

_cnt_a=$(samtools view -c "${BAM_A}"); echo "[FILT-A] ${PREFIX}: ${_cnt_a} reads → $(basename ${BAM_A})"
_cnt_b=$(samtools view -c "${BAM_B}"); echo "[FILT-B] ${PREFIX}: ${_cnt_b} reads → $(basename ${BAM_B})"
_cnt_c=$(samtools view -c "${BAM_C}"); echo "[FILT-C] ${PREFIX}: ${_cnt_c} reads → $(basename ${BAM_C})"
if [[ "${ALIGNER}" == "bwamem2" ]]; then
    echo "[HD-q20] ${PREFIX}: q20.fullRmsk → hd=0: $(samtools view -c -d hd:0 ${BAM_A})  hd=1: $(samtools view -c -d hd:1 ${BAM_A})"
fi

# ─── PER-SAMPLE ROW FILE ─────────────────────────────────────────────────────
_input_reads=$(samtools view -c "${BC_BAM}")
printf "%s\t%s\t%s\t%d\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%d\t%d\n" \
    "${PREFIX}" "${ALIGNER}" "${GENOME}" \
    "${_input_reads}" "${_trim_out}" "${_mapped}" "${_map_pct}" \
    "${_hd0}" "${_hd1}" "${_pct_hd1}" \
    "${_cnt_a}" "${_cnt_b}" "${_cnt_c}" \
    > "logs/mapUG_row_${PREFIX}.tsv"

echo ""
echo "════════════════════════════════════════════════════"
echo " 2.mapUG.sh — COMPLETE  $(date '+%Y-%m-%d %H:%M:%S')"
printf "  Sample   : %s\n" "${PREFIX}"
printf "  Raw BAM  : %s\n" "${RAW_BAM}"
printf "  Primary  : %s  [→ 3.run_bamsplitUG.sh]\n" "${BAM_A}"
echo "  Row file : logs/mapUG_row_${PREFIX}.tsv"
echo ""
echo "  After all tasks finish, assemble full summary:"
echo "  { printf 'sample\taligner\tgenome\tinput_reads(BCtagged)\ttrimmed_reads\tmapped_reads\tmapping_pct\thd0_mapped\thd1_mapped\thd1_pct\tq20_fullRmsk\tq0_fullRmsk\tq0_noRmsk\n'; cat logs/mapUG_row_*.tsv | sort; } > logs/mapUG_summary.tsv"
echo "════════════════════════════════════════════════════"
