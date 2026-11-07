#!/bin/bash
#SBATCH --job-name=AMBER_BamSplitUG
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --array=1-9
#SBATCH --output=BamSplit_slurmlogs/slurm-%A_%a.out

#Date: 2026-05-01
#Purpose: AMBER BamSplitting pipeline — Ultima UG edition.
#         Filtered BAM → UMI-tools group → BX sort → span-dedup + partition → cluster text + BAMs + summary.
#         CT mode ON by default (Cis/Trans/CT per size bin), BAM output ON by default.
#         NOTE: No PE flag removal [vs Illumina run_bamsplit.sh step 1/4] — UG reads are SE.
#Usage:
#  cd /path/to/data   (or Barnyard/)
#  sbatch --array=1-$(ls *.bwamem2.q20.fullRmsk.bam | wc -l) /path/to/AMBER/Processing/2.BamSplit/3.run_bamsplitUG.sh
#  sbatch --array=1-$(ls *.bwamem2.q0.fullRmsk.bam  | wc -l) /path/to/AMBER/Processing/2.BamSplit/3.run_bamsplitUG.sh "" q0_fullRmsk
#  sbatch --array=1-$(ls *.bwamem2.q0.noRmsk.bam    | wc -l) /path/to/AMBER/Processing/2.BamSplit/3.run_bamsplitUG.sh "" q0_noRmsk
#  Arg $1 (optional): bin spec string. Default: "1,2,3-5,6-10,11-30,31+"
#  Arg $2 (optional): filter level. One of: q20_fullRmsk (default) | q0_fullRmsk | q0_noRmsk
#  Bin format: exact int (1), range (3-5), or open-ended (31+). Comma-separated. No spaces inside.
#  Aligner is hardcoded to bwamem2 in CONFIG (toggle there if needed for bt2 comparison runs).
#  Run from the directory containing the filtered BAMs (e.g. Cortex/Process/).
#INPUT:  *.bwamem2.{filter}.bam   (from 2.mapUG.sh); filter level selected via $2
#OUTPUT: BamSplit/{FILTER_NAME}/{PREFIX}.{bin}_{CT}.bam + *_clusters.txt + final_summary.txt
#REQUIRES: module load samtools umi-tools (on login node before sbatch); python3 + pysam (pip install pysam)

set -eo pipefail

# ─── MODULES (load on login node before sbatch) ───────────────────────────────
module load samtools 2>/dev/null || true
module load umi-tools 2>/dev/null || true

# ─── CONFIG ───────────────────────────────────────────────────────────────────
BINS="${1:-1,2,3-5,6-10,11-30,31+}"   # override via CLI arg: sbatch ... "1-10,11-50,51+"
ALIGNER="bwamem2"                      # hardcoded; change here for bt2 comparison runs
THREADS=16

# Filter level set by $2 (default: q20_fullRmsk)
case "${2:-q20_fullRmsk}" in
  q0_fullRmsk) FILTER_NAME="q0_fullRmsk"; INPUT_SUFFIX=".${ALIGNER}.q0.fullRmsk.bam" ;;
  q0_noRmsk)   FILTER_NAME="q0_noRmsk";   INPUT_SUFFIX=".${ALIGNER}.q0.noRmsk.bam"   ;;
  *)           FILTER_NAME="q20_fullRmsk"; INPUT_SUFFIX=".${ALIGNER}.q20.fullRmsk.bam" ;;
esac

SPAN_SIZE=147       # span-dedup distance. 147 = nucleosome footprint (default), 200 = DNA-AMBER alt
EDIT_DIST=2         # UMI-tools hamming distance threshold (manuscript default: 2)
CT_FLAG="--ct"      # Cis/Trans splitting ON by default (clear to "" to disable)
RNA_FLAG=""         # set to "--rna" for RNA-AMBER merged BAMs
WRITE_BAMS="--write-bams"   # BAM output ON by default (clear to "" to skip BAMs)
# ─────────────────────────────────────────────────────────────────────────────

# Default: next to this script. Under sbatch the script is copied to a spool dir, so if you submit with sbatch
# run:  export BAMSPLIT_PY=/path/to/AMBER/Processing/2.BamSplit/2.amber_bamsplitUG.py   before sbatch.
BAMSPLIT_PY="${BAMSPLIT_PY:-$(dirname "$(readlink -f "$0")")/2.amber_bamsplitUG.py}"
[[ -f "${BAMSPLIT_PY}" ]] || { echo "ERROR: bamsplit script not found at ${BAMSPLIT_PY}. export BAMSPLIT_PY=/path/to/2.amber_bamsplitUG.py"; exit 1; }
OUTDIR="BamSplit/${FILTER_NAME}"

mkdir -p BamSplit_slurmlogs "${OUTDIR}"
[[ -n "${SLURM_ARRAY_TASK_ID}" ]] || { echo "ERROR: must be submitted via sbatch --array=..."; exit 1; }
INPUT_FILES=(*${INPUT_SUFFIX})
INPUT_BAM=${INPUT_FILES[$SLURM_ARRAY_TASK_ID-1]}
[[ -f "${INPUT_BAM}" ]] || { echo "ERROR: no BAM found for task ${SLURM_ARRAY_TASK_ID} (looked for ${INPUT_SUFFIX})"; exit 1; }
PREFIX=${INPUT_BAM%${INPUT_SUFFIX}}

UMI_BAM="${OUTDIR}/${PREFIX}.grouped.unsorted.bam"
SORTED_BAM="${OUTDIR}/${PREFIX}.grouped.sorted.bam"
OUTPUT_PREFIX="${OUTDIR}/${PREFIX}"

echo "════════════════════════════════════════════════════"
echo " 3.run_bamsplitUG.sh  $(date '+%Y-%m-%d %H:%M:%S')"
printf "  Job     : %-12s  Node: %s\n" "${SLURM_JOB_ID:-local}" "$(hostname -s)"
printf "  Aligner : %s   Span: %dbp   Edit-dist: %d   CT: %s\n" \
    "${ALIGNER}" "${SPAN_SIZE}" "${EDIT_DIST}" "${CT_FLAG:-(off)}"
printf "  Bins    : %s\n" "${BINS}"
printf "  Filter  : %s → %s\n" "${FILTER_NAME}" "${OUTDIR}"
printf "  Input   : %s\n" "${INPUT_BAM}"
echo "════════════════════════════════════════════════════"
echo ""

# [1/3] UMI-tools group (directional hamming-distance collapsing on 15bp SpatBC)
echo "[1/3] UMI-tools group (edit-distance=${EDIT_DIST})..."
umi_tools group -I "${INPUT_BAM}" --extract-umi-method=read_id --umi-separator='_' \
    --output-bam -S "${UMI_BAM}" --group-out "${OUTDIR}/${PREFIX}.groups.tsv" \
    -L "${OUTDIR}/${PREFIX}.umi.log" \
    --method=directional --edit-distance-threshold "${EDIT_DIST}"

# [2/3] Sort by BX tag; delete unsorted BAM immediately to save space
echo "[2/3] Sorting by BX tag..."
samtools sort -@ "${THREADS}" -t BX -o "${SORTED_BAM}" "${UMI_BAM}"
rm -f "${UMI_BAM}"

# [3/3] Span-dedup + partition → {bin}_{CT}.bam + {bin}_{CT}_clusters.txt + summary
python3 -c "import pysam" 2>/dev/null || { echo "ERROR: pysam not importable by $(which python3). Activate conda env with pysam before sbatch."; exit 1; }
echo "[3/3] Partitioning and generating clusters..."
python3 "${BAMSPLIT_PY}" -i "${SORTED_BAM}" -o "${OUTPUT_PREFIX}" \
    --bins "${BINS}" --span-size "${SPAN_SIZE}" -t "${THREADS}" ${CT_FLAG} ${RNA_FLAG} ${WRITE_BAMS} \
    | tee "${OUTPUT_PREFIX}.final_summary.txt"

# Delete large intermediates — not needed after bamsplit.
# To recreate: umi_tools group ... -S grouped.unsorted.bam → samtools sort -t BX → grouped.sorted.bam
rm -f "${SORTED_BAM}" "${OUTDIR}/${PREFIX}.groups.tsv"

echo ""
echo "════════════════════════════════════════════════════"
echo " 3.run_bamsplitUG.sh — COMPLETE  $(date '+%Y-%m-%d %H:%M:%S')"
echo "  Output dir   : ${OUTDIR}/"
echo "  Cluster files: ${OUTPUT_PREFIX}.{bin}[_{CT}]_clusters.txt"
echo "  BAM files    : ${OUTPUT_PREFIX}.{bin}[_{CT}].bam"
echo "  Summary      : ${OUTPUT_PREFIX}.final_summary.txt"
echo "  UMI log      : ${OUTDIR}/${PREFIX}.umi.log"
echo "════════════════════════════════════════════════════"
