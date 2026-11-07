#!/usr/bin/env python3
#Date: 2026-05-01
#Purpose: AMBER BamSplitting — Ultima UG edition. Single-pass through BX-sorted BAM:
#         span-dedup → partition by user-defined size bins (+ optional Cis/Trans) → cluster text/BAMs + summary.
#         Size bins fully configurable via --bins (no hardcoded singlet/pair/Clst labels).
#         hd:i clamp-quality tracked per-read in summary (bwamem2 only; absent in bt2).
#Usage:
#  python3 2.amber_bamsplitUG.py -i sample.grouped.sorted.bam -o outdir/sample
#  python3 2.amber_bamsplitUG.py -i ... --bins "1,2,3-5,6-10,11-30,31+" --ct --write-bams
#  python3 2.amber_bamsplitUG.py -i ... --bins "1-10,11-50,51-100,101+"   # custom bins
#  Bin format: exact int (1), range (3-5), or open-ended (31+). Comma-separated.
#INPUT:  BAM sorted by BX tag (output of: umi_tools group → samtools sort -t BX)
#OUTPUT: {prefix}.{bin}[_{Cis|Trans|CT}].bam  (if --write-bams)
#        {prefix}.{bin}[_{Cis|Trans|CT}]_clusters.txt
#        {prefix}.All[_{Cis|Trans|CT}].bam + _clusters.txt
#        Summary printed to stdout (tee'd to {prefix}.final_summary.txt by run_bamsplitUG.sh)

import argparse, sys
try:
    import pysam
except ImportError:
    sys.exit("Error: pysam required (pip install pysam)")

def parse_bins(bins_str):
    """Parse '1,2,3-5,6-10,11-30,31+' → list of (label, lo, hi) where hi=None means unbounded."""
    bins = []
    for tok in bins_str.split(','):
        tok = tok.strip()
        if tok.endswith('+'):
            bins.append((tok, int(tok[:-1]), None))
        elif '-' in tok:
            lo, hi = tok.split('-', 1)
            bins.append((tok, int(lo), int(hi)))
        else:
            n = int(tok); bins.append((tok, n, n))
    return bins

def get_size_bin(n, bins):
    for label, lo, hi in bins:
        if (hi is None and n >= lo) or (hi is not None and lo <= n <= hi): return label
    return bins[-1][0]   # fallback: assign to last bin if somehow unmatched

def get_hd(read):
    try: return read.get_tag('hd')
    except KeyError: return -1

def span_dedup(reads, span_size):
    reads.sort(key=lambda r: (r.reference_name, r.reference_start))
    last, out = {}, []
    for r in reads:
        c, s = r.reference_name, r.reference_start
        if c in last and (s - last[c]) <= span_size: continue
        last[c] = s; out.append(r)
    return out

def format_locus(r, rna):
    strand = '-' if r.flag & 16 else '+'
    chrom = r.reference_name
    if rna and chrom.startswith('RNA_'):
        return f"RNA[{strand}]_{chrom[4:]}:{r.reference_start+1}-{r.reference_end}"
    return f"DNA[{strand}]_{chrom}:{r.reference_start+1}-{r.reference_end}"

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-i", "--input-bam", required=True)
    p.add_argument("-o", "--output-prefix", required=True)
    p.add_argument("--bins", default="1,2,3-5,6-10,11-30,31+",
                   help="Size bin spec, e.g. '1,2,3-5,6-10,11-30,31+' or '1-10,11-50,51+'. "
                        "Each token: exact int (1), range (3-5), or open-ended (31+).")
    p.add_argument("--ct", action="store_true", help="Cis/Trans splitting (adds _Cis/_Trans/_CT per bin)")
    p.add_argument("--rna", action="store_true", help="RNA-AMBER mode (RNA_ prefix in chrom names)")
    p.add_argument("--write-bams", action="store_true", help="Output partitioned BAM files per bin")
    p.add_argument("--span-size", type=int, default=147, help="Span-dedup distance in bp (default: 147 = nucleosome footprint)")
    p.add_argument("-t", "--threads", type=int, default=8)
    args = p.parse_args()

    bins = parse_bins(args.bins)
    bin_labels = [b[0] for b in bins]
    print(f"--- amber_bamsplitUG.py | bins={args.bins} | span={args.span_size}bp | ct={'ON' if args.ct else 'OFF'} ---", file=sys.stderr)

    if args.ct:
        types = ["Cis", "Trans", "CT"]
        out_keys = [f"{lb}_{t}" for lb in bin_labels for t in types] + [f"All_{t}" for t in types]
    else:
        out_keys = bin_labels + ["All"]

    txt_handles = {k: open(f"{args.output_prefix}.{k}_clusters.txt", 'w') for k in out_keys}
    bam_handles = {}
    counts = {k: {"clusters": 0, "reads": 0} for k in out_keys}
    hd_total = {0: 0, 1: 0, -1: 0}
    n_raw_reads, n_raw_clusters = 0, 0

    try:
        with pysam.AlignmentFile(args.input_bam, "rb", threads=args.threads) as bam_in:
            header = bam_in.header.to_dict()
            if args.write_bams:
                bam_handles = {k: pysam.AlignmentFile(f"{args.output_prefix}.{k}.bam", "wb",
                               header=header, threads=args.threads) for k in out_keys}

            current_bx, cluster = None, []

            def flush_cluster():
                nonlocal n_raw_reads, n_raw_clusters
                if not cluster: return
                n_raw_reads += len(cluster); n_raw_clusters += 1
                deduped = span_dedup(cluster, args.span_size)
                if not deduped: return
                for r in deduped:
                    hd = get_hd(r); hd_total[hd if hd in (0, 1) else -1] += 1
                sb = get_size_bin(len(deduped), bins)
                bx = cluster[0].get_tag('BX')
                loci = [format_locus(r, args.rna) for r in deduped]
                line = f"{bx} {' '.join(loci)}\n"
                n = len(deduped)
                if args.ct:
                    chroms = {(r.reference_name[4:] if args.rna and r.reference_name.startswith('RNA_')
                               else r.reference_name) for r in deduped}
                    ct = "Cis" if len(chroms) == 1 else "Trans"
                    dest_keys = [f"{sb}_{ct}", f"{sb}_CT", f"All_{ct}", "All_CT"]
                else:
                    dest_keys = [sb, "All"]
                for k in dest_keys:
                    txt_handles[k].write(line); counts[k]["clusters"] += 1; counts[k]["reads"] += n
                    if k in bam_handles:
                        for r in deduped: bam_handles[k].write(r)

            for read in bam_in:
                try: bx = read.get_tag('BX')
                except KeyError: continue
                if bx != current_bx:
                    flush_cluster(); cluster, current_bx = [read], bx
                else: cluster.append(read)
            flush_cluster()

    finally:
        for h in txt_handles.values(): h.close()
        for h in bam_handles.values(): h.close()

    pct = lambda x, t: f"({100*x/t:.2f}%)" if t > 0 else "(0.00%)"

    if args.ct:
        total_c = counts["All_Cis"]["clusters"] + counts["All_Trans"]["clusters"]
        total_r = counts["All_Cis"]["reads"] + counts["All_Trans"]["reads"]
        print(f"\n{'='*72}\nAMBER BamSplitting Summary (UG) | span={args.span_size}bp | bins={args.bins} | CT mode\n{'='*72}")
        print(f"Total Reads (Unique Loci): {total_r:,}\nTotal Final Clusters: {total_c:,}\n")
        hdr = f"{'Bin':<12} {'Type':<6} {'Clusters':>12} {'(%)':>9}  {'Reads':>12} {'(%)':>9}"
        print(hdr); print("-" * len(hdr))
        for lb in bin_labels:
            for t in ["Cis", "Trans", "CT"]:
                k = f"{lb}_{t}"; c, r = counts[k]["clusters"], counts[k]["reads"]
                print(f"{lb:<12} {t:<6} {c:>12,} {pct(c, total_c):>9}  {r:>12,} {pct(r, total_r):>9}")
            print(f"{'':->68}")
    else:
        total_c = counts["All"]["clusters"]; total_r = counts["All"]["reads"]
        print(f"\n{'='*68}\nAMBER BamSplitting Summary (UG) | span={args.span_size}bp | bins={args.bins}\n{'='*68}")
        print(f"Total Reads (Unique Loci): {total_r:,}\nTotal Final Clusters: {total_c:,}\n")
        for lb in bin_labels:
            c, r = counts[lb]["clusters"], counts[lb]["reads"]
            print(f"{lb:<18} {c:>12,} {pct(c, total_c):>9}    {r:>12,} {pct(r, total_r):>9}")

    hd_known = hd_total[0] + hd_total[1]
    if hd_known > 0:
        print(f"\n── hd:i clamp quality (post-dedup reads) ──────────────────────────")
        print(f"  hd=0 (exact clamp):      {hd_total[0]:>10,}  {pct(hd_total[0], hd_known)}")
        print(f"  hd=1 (1-mismatch clamp): {hd_total[1]:>10,}  {pct(hd_total[1], hd_known)}")
    elif hd_total[-1] > 0:
        print(f"\n── hd:i clamp quality: n/a (bt2 mode — tag not propagated) ────────")

    n_removed = n_raw_reads - total_r; n_dropped_clusters = n_raw_clusters - total_c
    print(f"\n── Span-dedup efficiency ({args.span_size}bp window) ──────────────────────────────")
    print(f"  Input reads (post-UMI-group):   {n_raw_reads:>13,}")
    print(f"  Reads retained (unique loci):   {total_r:>13,}  ({100*total_r/n_raw_reads:.1f}% of input)" if n_raw_reads > 0 else f"  Reads retained (unique loci):   {total_r:>13,}")
    print(f"  Reads removed by span-dedup:    {n_removed:>13,}  ({100*n_removed/n_raw_reads:.1f}% removed)" if n_raw_reads > 0 else f"  Reads removed by span-dedup:    {n_removed:>13,}")
    print(f"  Input BX groups:                {n_raw_clusters:>13,}")
    print(f"  BX groups retained:             {total_c:>13,}  ({n_dropped_clusters:,} dropped — all loci within span)")

    n_files = sum(1 for k in out_keys if counts[k]["clusters"] > 0)
    print(f"\nWrote {n_files} non-empty cluster files → {args.output_prefix}.{{bin}}_clusters.txt")
    if args.write_bams: print(f"Wrote {n_files} partitioned BAMs → {args.output_prefix}.{{bin}}.bam")

if __name__ == "__main__":
    main()
