#!/bin/bash
#Date: 2026-05-26
#Purpose: Apply CT size cutoff to a cluster file.
#         Single cutoff (CT >= N) or range (CT >= LO and CT <= HI).
#         CT = number of loci per SpatBC = NF-1 (line has SpatBC + N loci).
#Usage: bash /path/to/AMBER/Processing/2.BamSplit/4.SizeCutoff.sh \
#         <*_clusters.txt> <cutoff_lo> [cutoff_hi]
#       Single:  bash 4.SizeCutoff.sh Master_All_CT_clusters.txt 11        -> CT11co.txt
#       Range:   bash 4.SizeCutoff.sh Master_All_CT_clusters.txt 3 10      -> CT3-10co.txt
#       Multiple single: for co in 11 15; do bash 4.SizeCutoff.sh Master_All_CT_clusters.txt $co; done
#INPUT: *_clusters.txt from 2.amber_bamsplitUG.py
#       Line format: "SPATBC locus1 locus2 ..." — CT = NF-1 (fields minus the SpatBC column)

set -uo pipefail

[[ $# -lt 2 ]] && { echo "Usage: $0 <clusters.txt> <cutoff_lo> [cutoff_hi]"; exit 1; }

IN="$1"
LO="$2"
HI="${3:-}"

if [[ -n "$HI" ]]; then
    OUT="${IN%.txt}_CT${LO}-${HI}co.txt"
    awk -v lo="$LO" -v hi="$HI" 'NF-1 >= lo && NF-1 <= hi' "$IN" > "$OUT"
    LABEL="CT ${LO}–${HI}"
else
    OUT="${IN%.txt}_CT${LO}co.txt"
    awk -v lo="$LO" 'NF-1 >= lo' "$IN" > "$OUT"
    LABEL="CT≥${LO}"
fi

N_IN=$(wc -l < "$IN")
N_OUT=$(wc -l < "$OUT")
echo "Input:  $IN ($N_IN lines)"
echo "Output: $OUT ($N_OUT lines, ${LABEL})"
awk "BEGIN{printf \"Retained: %.4f%%\n\", $N_OUT/$N_IN*100}"
