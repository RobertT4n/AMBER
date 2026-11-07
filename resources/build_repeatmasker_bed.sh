#!/bin/bash
#Date: 2026-11-07
#Purpose: Build the RepeatMasker bed used by 2.mapUG.sh to remove reads in repeats (UCSC RepeatMasker, milliDiv < 140).
#Usage:   bash build_repeatmasker_bed.sh <hg38|mm10> [OUTPUT.bed]
#         Combined human+mouse bed (for the hgmm species-mixing genome):
#           bash build_repeatmasker_bed.sh hg38 hg38.bed; bash build_repeatmasker_bed.sh mm10 mm10.bed
#           { awk 'BEGIN{OFS="\t"}{$1="hg_"$1; print}' hg38.bed; awk 'BEGIN{OFS="\t"}{$1="mm_"$1; print}' mm10.bed; } > hg38_mm10_combined.bed
#INPUT:   UCSC rmsk table (downloaded): https://hgdownload.soe.ucsc.edu/goldenPath/<assembly>/database/rmsk.txt.gz
#OUTPUT:  6-column bed: chrom  start  end  repName  milliDiv  strand   (sorted)
#REQUIRES: curl (or wget), gzip, awk, sort

set -euo pipefail
[[ $# -ge 1 ]] || { echo "Usage: $0 <hg38|mm10> [OUTPUT.bed]"; exit 1; }
ASM="$1"
case "${ASM}" in hg38|mm10) ;; *) echo "ERROR: assembly must be hg38 or mm10"; exit 1;; esac
OUT="${2:-${ASM}.repeatmasker.milliDiv_lt140.bed}"
[[ ! -e "${OUT}" ]] || { echo "ERROR: ${OUT} already exists; not overwriting"; exit 1; }

URL="https://hgdownload.soe.ucsc.edu/goldenPath/${ASM}/database/rmsk.txt.gz"
echo "Downloading ${URL}"
# rmsk.txt columns (1-based): 1 bin, 2 swapScore, 3 milliDiv, 4 milliDel, 5 milliIns, 6 genoName, 7 genoStart, 8 genoEnd,
#                             9 genoLeft, 10 strand, 11 repName, 12 repClass, 13 repFamily, ...
curl -sSL "${URL}" | gzip -dc | \
  awk -F'\t' 'BEGIN{OFS="\t"} $3 < 140 {print $6, $7, $8, $11, $3, $10}' | \
  sort -k1,1 -k2,2n > "${OUT}"

echo "Wrote ${OUT}: $(wc -l < "${OUT}") repeats with milliDiv < 140"

# ---- Provenance ----
# Upstream:   UCSC Genome Browser rmsk table (RepeatMasker annotations), see URL above
# Downstream: 1.Raw2Bam/2.mapUG.sh (bedtools intersect -v removes reads overlapping these intervals)
# References: UCSC RepeatMasker track; milliDiv < 140 repeat-masking threshold
# Notes:      milliDiv = per-mille divergence from the repeat consensus. Only repeats with milliDiv < 140 are written (the mask);
#             the most diverged (oldest) repeats, milliDiv >= 140, are not masked.
#             Not yet cross-checked line-by-line against the bed originally used for the manuscript (compare line counts before relying on it).
