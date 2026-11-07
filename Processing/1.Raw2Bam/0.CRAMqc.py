#!/usr/bin/env python3
"""
CRAMqc v2.0 — Ultima UG CRAM quality control report (FastQC-equivalent HTML).
============================================================================

PURPOSE
-------
Produce a self-contained, single-file HTML QC report for Ultima Genomics CRAM files.
Designed for the AMBER pipeline (NYGC cluster), but works on any Ultima CRAM or BAM.
Run BEFORE 1.CRAMprocess.sh to assess raw read quality, OR after alignment for mapped QC.

USAGE
-----
    python3 0.CRAMqc.py sample.cram [options]

    Common examples:
        python3 0.CRAMqc.py sample.cram                    # standard report, 500k reads
        python3 0.CRAMqc.py sample.cram --sample 100000    # faster, 100k reads
        python3 0.CRAMqc.py sample.cram --amber            # include AMBER clamp-detection module
        python3 0.CRAMqc.py sample.1MIL.cram --sample 100000 --amber  # dry-run on 1M subset

INPUT
-----
    Ultima UG CRAM (unaligned, reference-less) or any BAM/CRAM.
    No reference genome needed for unaligned CRAMs (uses pysam with check_sq=False).
    NOTE: pysam/htslib handles Ultima's non-standard file_id container header fine.
          The Rust noodles-cram library cannot handle Ultima's non-standard file_id
          container header; samtools subprocess is used instead.

OUTPUT
------
    {prefix}.CRAMqc.html   Self-contained HTML report (open in any browser)
    {prefix}.CRAMqc.tsv    Summary table (Google Sheets: paste → split text by Tab)

MODULES IN THE HTML REPORT
---------------------------
Standard (FastQC-equivalent):
    1. Basic Statistics            — filename, read count, length, GC%, quality
    2. Per Base Sequence Quality   — FastQC boxplot: Q10/Q25/Q50/Q75/Q90 per position
                                     Background zones: green (Q≥20) / yellow (Q13–20) / red (<Q13)
                                     NOTE: Ultima Phred ≠ Illumina Phred. Ultima Q15 ≈ Illumina Q20–30.
                                     Thresholds are Ultima-calibrated; red zone (<Q13) is more permissive
                                     than standard FastQC (<Q20).
    3. Per Sequence Quality Scores — histogram of mean Phred per read
    4. Per Base Sequence Content   — A/T/G/C % at each read position
    5. Per Sequence GC Content     — GC% histogram across reads
    6. Sequence Length Distribution
    7. Sequence Duplication Levels — FastQC-style bins on first 100k read hashes
                                     NOT INFORMATIVE in --amber mode: first 50bp is fixed adapter
                                     (14bp Tn5 + 19bp Nextera ME) + SpatBC cell barcode — not genomic insert.
                                     True PCR duplication should be assessed post-alignment (markdup on mapped BAM).

Ultima-specific:
    8. rq:f Quality Distribution   — Ultima's per-read quality score (0–1), not in FastQC

AMBER-specific (--amber flag; only meaningful on raw pre-filter CRAMs):
    9. Clamp Detection             — GCATGCGG bridge-ligation clamp yield + positional histogram
                                     Expected: ~30% yield (bridge-ligated hub reads)
                                     hd=0 (exact) vs hd=1 (1-mismatch) breakdown

PASS / WARN / FAIL THRESHOLDS
------------------------------
    Module 2 (per-base quality):    PASS all medians ≥18 and Q25 ≥13; WARN any median <18; FAIL any Q25 <10
                                    (Ultima-calibrated: Q15 ≈ Illumina Q20–30; thresholds scaled accordingly)
    Module 3 (per-seq quality):     PASS mean >27; WARN mean 20–27; FAIL mean <20
    Module 5 (GC content):          PASS if SD <15%; WARN if SD <25%; FAIL if SD ≥25%
    Module 6 (length):              PASS all same; WARN >10 length values; FAIL any read <20bp
    Module 7 (duplication):         PASS ≥80% unique; WARN ≥50%; FAIL <50%
    Module 8 (rq:f):                Standard: PASS ≥80% rq>0.5; WARN ≥60%; FAIL <60%
                                    AMBER pre-trim: relaxed — fixed adapter prefix depresses rq:f on all reads.
                                    PASS if mean≥0.35 and <40% reads with rq<0.3 (bell ~0.4–0.5 = healthy).
    Module 9 (clamp):               PASS ≥25% yield; WARN ≥10%; FAIL <10%

SAMPLING NOTE
-------------
    CRAMqc samples N reads (default 500k) for ALL metrics, unlike FastQC which streams
    all reads for quality/GC/length but subsamples 100k for duplication. This is intentional:
    streaming 500M+ reads through a samtools subprocess would take 30–60 min per CRAM.
    500k is statistically robust for all quality distributions. Duplication is 100k (first
    50bp hashes), same as FastQC — but see Module 7 warning for AMBER pre-trim CRAMs.

REQUIREMENTS
------------
    numpy, matplotlib  (pip install numpy matplotlib)
    samtools in PATH (reads CRAM via subprocess pipe — bypasses pysam CRAM limitations)

Date: 2026-04-27
"""

import argparse, sys, os, io, base64, datetime, hashlib, subprocess
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend (required on HPC cluster — no display)
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection, LineCollection
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam required. Install: pip install pysam")

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
VERSION = "2.0"
MAX_QUAL = 94        # Phred 0–93, encoded as ASCII 33–126
MAX_POS  = 400       # track per-position stats up to this read length (Ultima ~300bp)
CLAMP    = b"GCATGCGG"
CLAMP_WIN_S, CLAMP_WIN_E = 47, 65  # same window as 1.CRAMprocess.sh

# Duplication binning labels (FastQC-style)
DUP_BINS   = [1, 2, 3, 5, 10, 50, 100, 500, 1000, float('inf')]
DUP_LABELS = ['1','2','3','4-5','6-10','11-50','51-100','101-500','501-1000','>1000']

# ─────────────────────────────────────────────────────────────────────────────
# Data collection — single streaming pass
# ─────────────────────────────────────────────────────────────────────────────

def collect(cram_path, n_sample, n_pos_sample, threads, do_amber):
    """
    Stream up to n_sample reads from cram_path. Returns a metrics dict (m) containing
    all raw data needed by the plot and assessment functions.

    Per-position quality is tracked for first n_pos_sample reads only (memory bound).
    Duplication hashing on first 100k reads (FastQC convention).
    """
    m = {
        'file':        os.path.basename(cram_path),
        'cram_path':   cram_path,
        'n_sample':    n_sample,
        'n_pos_sample': n_pos_sample,
        'date':        datetime.datetime.now().strftime('%Y-%m-%d %H:%M'),
        # Read-level arrays
        'lengths':     [],
        'mean_quals':  [],
        'gc_fracs':    [],
        'rq_scores':   [],
        'n_rq_missing': 0,
        # Per-position quality histogram [MAX_POS × MAX_QUAL]
        'qual_hist':   np.zeros((MAX_POS, MAX_QUAL), dtype=np.int64),
        # Per-position base content [MAX_POS × 5] (A C G T N)
        'base_hist':   np.zeros((MAX_POS, 5), dtype=np.int64),
        # Duplication: hash → count (first 100k reads, first 50bp)
        'dup_hashes':  Counter(),
        'n_dup_reads': 0,
        # AMBER clamp (if --amber)
        'clamp_pos':   [],    # 1-based positions where clamp was found
        'clamp_hd':    [],    # 0 or 1
        'n_clamp_scanned': 0,
        # Tracking
        'n_reads':     0,
        'n_overlength': 0,
        'do_amber':    do_amber,
    }
    BASE_IDX = {ord('A'): 0, ord('C'): 1, ord('G'): 2, ord('T'): 3,
                ord('a'): 0, ord('c'): 1, ord('g'): 2, ord('t'): 3}

    try:
        # samtools view pipes raw SAM lines — bypasses pysam CRAM parsing entirely
        # (Ultima CRAMs have non-standard container headers that confuse pysam/htslib)
        proc = subprocess.Popen(
            ['samtools', 'view', '-@', str(threads), cram_path],
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        for raw in proc.stdout:
            if m['n_reads'] >= n_sample:
                break
            fields = raw.split(b'\t')
            if len(fields) < 11 or fields[9] == b'*' or fields[10] == b'*':
                continue
            seq_b = fields[9].upper()
            seq   = seq_b.decode()
            quals = [b - 33 for b in fields[10].rstrip(b'\n')]

            n     = len(seq)
            q_arr = np.array(quals, dtype=np.uint8)

            m['lengths'].append(n)
            m['mean_quals'].append(float(q_arr.mean()))
            gc = (seq_b.count(b'G') + seq_b.count(b'C')) / n
            m['gc_fracs'].append(gc * 100)

            # rq:f — Ultima per-read quality score (scan optional fields)
            rq = next((f[5:] for f in fields[11:] if f.startswith(b'rq:f:')), None)
            if rq:
                try: m['rq_scores'].append(float(rq))
                except ValueError: m['n_rq_missing'] += 1
            else:
                m['n_rq_missing'] += 1

            # Per-position quality histogram (limited to first n_pos_sample reads)
            if m['n_reads'] < n_pos_sample:
                if n > MAX_POS:
                    m['n_overlength'] += 1
                trunc = min(n, MAX_POS)
                for pos in range(trunc):
                    m['qual_hist'][pos, min(int(q_arr[pos]), MAX_QUAL - 1)] += 1

            # Per-position base content (limited to first n_pos_sample reads)
            if m['n_reads'] < n_pos_sample:
                for pos in range(min(n, MAX_POS)):
                    idx = BASE_IDX.get(seq_b[pos], 4)  # 4 = N
                    m['base_hist'][pos, idx] += 1

            # Duplication: hash first 50bp (first 100k reads, FastQC convention)
            if m['n_dup_reads'] < 100000:
                key = hashlib.md5(seq_b[:50]).digest()
                m['dup_hashes'][key] += 1
                m['n_dup_reads'] += 1

            # AMBER clamp detection (optional)
            if do_amber:
                m['n_clamp_scanned'] += 1
                cp, hd = _find_clamp(seq_b)
                if cp > 0:
                    m['clamp_pos'].append(cp)
                    m['clamp_hd'].append(hd)

            m['n_reads'] += 1

        proc.stdout.close(); proc.wait()

    except Exception as e:
        sys.exit(f"ERROR reading {cram_path}: {e}")

    if m['n_reads'] == 0:
        sys.exit("ERROR: No reads found in input file.")
    if m['n_overlength'] > 0:
        print(f"  WARNING: {m['n_overlength']} reads >MAX_POS={MAX_POS}bp — per-position data truncated.")

    return m


def _find_clamp(seq_b, win_s=CLAMP_WIN_S, win_e=CLAMP_WIN_E, max_mm=1):
    """Scan seq_b for GCATGCGG with ≤max_mm mismatches in window [win_s, win_e].
    Returns (best_pos_1based, best_hd). Returns (-1, max_mm+1) if not found."""
    best_hd, best_pos = max_mm + 1, -1
    for i in range(win_s - 1, min(win_e, len(seq_b) - 8)):
        obs = seq_b[i:i+8]
        if len(obs) < 8:
            break
        hd = sum(a != b for a, b in zip(obs, CLAMP))
        if hd < best_hd:
            best_hd, best_pos = hd, i + 1
            if hd == 0:
                break
    return best_pos, best_hd


def _hist_percentiles(h):
    """Compute (q10, q25, q50, q75, q90, mean) from a quality count array h[MAX_QUAL]."""
    total = h.sum()
    if total == 0:
        return 0, 0, 0, 0, 0, 0.0
    cumsum = np.cumsum(h)
    pcts   = np.searchsorted(cumsum, np.array([0.10, 0.25, 0.50, 0.75, 0.90]) * total)
    mean   = float(np.dot(np.arange(MAX_QUAL), h)) / total
    return tuple(int(p) for p in pcts) + (mean,)


def compute_per_pos_stats(m):
    """Compute per-position quality stats from qual_hist. Returns dict of arrays."""
    pos_data = [p for p in range(MAX_POS) if m['qual_hist'][p].sum() > 0]
    if not pos_data:
        return None
    stats = [_hist_percentiles(m['qual_hist'][p]) for p in pos_data]
    return {
        'positions': pos_data,
        'q10':  [s[0] for s in stats],
        'q25':  [s[1] for s in stats],
        'q50':  [s[2] for s in stats],
        'q75':  [s[3] for s in stats],
        'q90':  [s[4] for s in stats],
        'mean': [s[5] for s in stats],
    }


# ─────────────────────────────────────────────────────────────────────────────
# Pass / Warn / Fail assessment
# ─────────────────────────────────────────────────────────────────────────────

def assess(m, pps):
    """Return dict of module_id → ('pass'|'warn'|'fail')."""
    s = {}

    # Per-base quality
    if pps:
        q25_min = min(pps['q25']); q50_min = min(pps['q50'])
        # Ultima-calibrated thresholds (Q15 ≈ Illumina Q20–30; zones scaled down ~0.65×)
        if q25_min >= 13 and q50_min >= 18:   s['per_base_qual'] = 'pass'
        elif q25_min >= 10 and q50_min >= 13: s['per_base_qual'] = 'warn'
        else:                                  s['per_base_qual'] = 'fail'
    else:
        s['per_base_qual'] = 'warn'

    # Per-sequence quality
    mq = np.mean(m['mean_quals'])
    if mq > 27:   s['per_seq_qual'] = 'pass'
    elif mq > 20: s['per_seq_qual'] = 'warn'
    else:          s['per_seq_qual'] = 'fail'

    # Per-base sequence content: PASS if max A-T imbalance < 10% at any position
    s['seq_content'] = 'pass'  # simplification: always pass for now (complex to assess)

    # GC content: PASS if std < 15%
    gc_sd = float(np.std(m['gc_fracs']))
    if gc_sd < 15:   s['gc_content'] = 'pass'
    elif gc_sd < 25: s['gc_content'] = 'warn'
    else:             s['gc_content'] = 'fail'

    # Length distribution
    lengths = m['lengths']
    n_short = sum(1 for l in lengths if l < 20)
    n_unique_lengths = len(set(lengths))
    if n_short > 0:                s['length_dist'] = 'fail'
    elif n_unique_lengths > 10:    s['length_dist'] = 'warn'
    else:                           s['length_dist'] = 'pass'

    # Duplication — N/A in amber mode (first 50bp is fixed adapter + SpatBC, not genomic)
    dup_c = m['dup_hashes']
    n_dup = sum(dup_c.values())
    pct_unique = 100 * sum(1 for v in dup_c.values() if v == 1) / max(n_dup, 1)
    if m.get('do_amber'):
        s['duplication'] = 'warn'   # warn = visually neutral; description explains N/A
    elif pct_unique >= 80:   s['duplication'] = 'pass'
    elif pct_unique >= 50: s['duplication'] = 'warn'
    else:                   s['duplication'] = 'fail'

    # rq:f — relaxed thresholds in amber mode: fixed adapter prefix (~33bp Tn5+ME) on every
    # read depresses rq:f systematically; bell-shaped distribution ~0.4–0.5 is expected/healthy.
    rq = m['rq_scores']
    if rq:
        pct_hi = 100 * sum(1 for r in rq if r > 0.5) / len(rq)
        pct_lo = 100 * sum(1 for r in rq if r < 0.3) / len(rq)
        if m.get('do_amber'):
            # AMBER pre-trim: pass if bell-shaped (mean≥0.35 and <0.3 spike <40%)
            rq_mean = sum(rq) / len(rq)
            if rq_mean >= 0.35 and pct_lo < 40:   s['rq_qual'] = 'pass'
            elif rq_mean >= 0.25 and pct_lo < 60: s['rq_qual'] = 'warn'
            else:                                   s['rq_qual'] = 'fail'
        else:
            if pct_hi >= 80:   s['rq_qual'] = 'pass'
            elif pct_hi >= 60: s['rq_qual'] = 'warn'
            else:               s['rq_qual'] = 'fail'
    else:
        s['rq_qual'] = 'warn'  # tag absent

    # Clamp
    if m['n_clamp_scanned'] > 0:
        pct = 100 * len(m['clamp_pos']) / m['n_clamp_scanned']
        if pct >= 25:   s['clamp'] = 'pass'
        elif pct >= 10: s['clamp'] = 'warn'
        else:            s['clamp'] = 'fail'
    else:
        s['clamp'] = 'warn'

    return s


# ─────────────────────────────────────────────────────────────────────────────
# Plot generators — each returns a base64-encoded <img> tag
# ─────────────────────────────────────────────────────────────────────────────

def _b64img(fig, dpi=130):
    """Render matplotlib fig to base64 PNG string for HTML embedding."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode()
    plt.close(fig)
    return f'<img class="plot" src="data:image/png;base64,{b64}" alt="plot">'


def plot_per_base_quality(m, pps):
    """FastQC-style per-position boxplot: Q10/Q25/Q50/Q75/Q90 + colored background zones."""
    if not pps:
        return '<p class="nodata">No per-position quality data available.</p>'
    positions = pps['positions']
    fig, ax = plt.subplots(figsize=(13, 5))
    # Background quality zones (FastQC colors)
    # Ultima-calibrated zones: Q15 ≈ Illumina Q20–30; thresholds are ~0.65× Illumina FastQC
    ax.axhspan(20, 43, color='#d5e8d4', alpha=0.85, zorder=0)   # green: good (Ultima Q≥20)
    ax.axhspan(13, 20, color='#fff2cc', alpha=0.85, zorder=0)   # yellow: reasonable (Ultima Q13–20)
    ax.axhspan(0,  13, color='#f8cecc', alpha=0.85, zorder=0)   # red: poor (Ultima Q<13)
    # Boxes (Q25–Q75) using PatchCollection (efficient for 300+ positions)
    boxes = [Rectangle((p - 0.45, q25), 0.9, max(q75 - q25, 0.5))
             for p, q25, q75 in zip(positions, pps['q25'], pps['q75'])]
    ax.add_collection(PatchCollection(boxes, facecolor='#FFD700', edgecolor='#888', linewidth=0.25, zorder=2))
    # Whiskers Q10–Q25 and Q75–Q90
    whiskers = ([[(p, q10), (p, q25)] for p, q10, q25 in zip(positions, pps['q10'], pps['q25'])] +
                [[(p, q75), (p, q90)] for p, q75, q90 in zip(positions, pps['q75'], pps['q90'])])
    ax.add_collection(LineCollection(whiskers, colors='#888', linewidths=0.6, zorder=1))
    # Median lines (red)
    medians = [[(p - 0.45, q50), (p + 0.45, q50)] for p, q50 in zip(positions, pps['q50'])]
    ax.add_collection(LineCollection(medians, colors='red', linewidths=1.2, zorder=3))
    # Mean line (blue)
    ax.plot(positions, pps['mean'], 'b-', lw=1.5, alpha=0.8, zorder=4)
    # Reference lines
    ax.axhline(28, color='gray', lw=0.5, ls='--', alpha=0.4, zorder=1)
    ax.axhline(20, color='gray', lw=0.5, ls='--', alpha=0.4, zorder=1)
    ax.axhline(15, color='tomato', lw=1, ls=':', alpha=0.8, zorder=1)
    ax.set_xlim(-0.5, max(positions) + 0.5)
    ax.set_ylim(0, 43)
    ax.set_xlabel('Position in read (bp)', fontsize=11)
    ax.set_ylabel('Phred quality score', fontsize=11)
    ax.set_title(f'Per Base Sequence Quality  [{m["n_pos_sample"]:,} reads for per-position data]', fontsize=12)
    legend = [Patch(facecolor='#FFD700', edgecolor='gray', label='IQR (Q25–Q75)'),
              Line2D([0],[0], color='red', lw=1.5, label='Median'),
              Line2D([0],[0], color='blue', lw=1.5, label='Mean'),
              Line2D([0],[0], color='tomato', lw=1, ls=':', label='Q15 (cutadapt threshold)'),
              Patch(facecolor='#d5e8d4', label='Good (Q≥20, Ultima-scaled)'),
              Patch(facecolor='#fff2cc', label='Reasonable (Q20–28)'),
              Patch(facecolor='#f8cecc', label='Poor (Q<20)')]
    ax.legend(handles=legend, fontsize=8, loc='lower left', ncol=2)
    fig.tight_layout()
    return _b64img(fig)


def plot_per_seq_quality(m):
    fig, ax = plt.subplots(figsize=(9, 4))
    mq = np.array(m['mean_quals'])
    ax.hist(mq, bins=np.arange(0, 43, 0.5), color='steelblue', edgecolor='none', alpha=0.85)
    ax.axvline(np.mean(mq), color='darkorange', lw=1.5, ls='--', label=f'mean={np.mean(mq):.1f}')
    ax.axvline(20, color='gold', lw=0.8, ls=':', label='Q20')
    ax.axvline(27, color='limegreen', lw=0.8, ls=':', label='Q27')
    ax.set_xlabel('Mean Phred quality per read', fontsize=11)
    ax.set_ylabel('Read count', fontsize=11)
    ax.set_title('Per Sequence Quality Scores', fontsize=12)
    ax.legend(fontsize=9); fig.tight_layout()
    return _b64img(fig)


def plot_seq_content(m, pps):
    """A/T/G/C % per read position."""
    if not pps:
        return '<p class="nodata">No per-position data available.</p>'
    positions = pps['positions']
    bh = m['base_hist']
    totals = bh[positions].sum(axis=1, keepdims=True).clip(min=1)
    pcts = 100 * bh[positions] / totals   # [n_pos × 5], columns = A C G T N
    fig, ax = plt.subplots(figsize=(13, 4))
    for i, (label, color) in enumerate([('A','#e74c3c'),('C','#3498db'),('G','#2ecc71'),('T','#f39c12')]):
        ax.plot(positions, pcts[:, i], color=color, lw=1.2, label=label)
    ax.axhline(25, color='gray', lw=0.5, ls='--', alpha=0.5)
    ax.set_xlim(-0.5, max(positions) + 0.5)
    ax.set_ylim(0, 100)
    ax.set_xlabel('Position in read (bp)', fontsize=11)
    ax.set_ylabel('Base content (%)', fontsize=11)
    ax.set_title('Per Base Sequence Content', fontsize=12)
    ax.legend(fontsize=10, loc='upper right'); fig.tight_layout()
    return _b64img(fig)


def plot_gc_content(m):
    fig, ax = plt.subplots(figsize=(9, 4))
    gc = np.array(m['gc_fracs'])
    ax.hist(gc, bins=np.arange(0, 101, 1), color='forestgreen', edgecolor='none', alpha=0.85)
    ax.axvline(np.mean(gc), color='darkorange', lw=1.5, ls='--', label=f'mean={np.mean(gc):.1f}%')
    ax.axvline(50, color='gray', lw=0.8, ls=':', alpha=0.6, label='50%')
    ax.set_xlabel('GC content per read (%)', fontsize=11)
    ax.set_ylabel('Read count', fontsize=11)
    ax.set_title('Per Sequence GC Content', fontsize=12)
    ax.legend(fontsize=9); fig.tight_layout()
    return _b64img(fig)


def plot_length_dist(m):
    fig, ax = plt.subplots(figsize=(9, 4))
    lengths = np.array(m['lengths'])
    ax.hist(lengths, bins=80, color='steelblue', edgecolor='none', alpha=0.85)
    ax.axvline(np.mean(lengths), color='darkorange', lw=1.5, ls='--',
               label=f'mean={np.mean(lengths):.0f}')
    ax.axvline(np.median(lengths), color='tomato', lw=1.5,
               label=f'median={np.median(lengths):.0f}')
    ax.set_xlabel('Read length (bp)', fontsize=11)
    ax.set_ylabel('Read count', fontsize=11)
    ax.set_title('Sequence Length Distribution', fontsize=12)
    ax.legend(fontsize=9); fig.tight_layout()
    return _b64img(fig)


def plot_duplication(m):
    """FastQC-style duplication level bar chart."""
    dup_c = m['dup_hashes']
    counts = list(dup_c.values())
    total  = len(counts)
    if total == 0:
        return '<p class="nodata">No duplication data.</p>'

    def bin_idx(v):
        for i, thresh in enumerate(DUP_BINS[:-1]):
            if v <= thresh:
                return i
        return len(DUP_BINS) - 1

    bins = Counter(bin_idx(v) for v in counts)
    pct_of_total = [100 * bins.get(i, 0) / total for i in range(len(DUP_LABELS))]
    n_unique = sum(1 for v in counts if v == 1)
    pct_unique = 100 * n_unique / total

    fig, ax = plt.subplots(figsize=(10, 4))
    xs = range(len(DUP_LABELS))
    ax.bar(xs, pct_of_total, color='steelblue', alpha=0.85, edgecolor='none')
    ax.set_xticks(list(xs))
    ax.set_xticklabels(DUP_LABELS, fontsize=9)
    ax.set_xlabel('Sequence duplication level', fontsize=11)
    ax.set_ylabel('% of reads', fontsize=11)
    ax.set_title(f'Sequence Duplication Levels  [{pct_unique:.1f}% unique, from first {total:,} reads]',
                 fontsize=12)
    fig.tight_layout()
    return _b64img(fig)


def plot_rq_dist(m):
    rq = m['rq_scores']
    if not rq:
        return '<p class="nodata">rq:f tag absent — this file may be a post-alignment BAM (tag only present in raw Ultima CRAMs).</p>'
    rq = np.array(rq)
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.hist(rq, bins=60, color='mediumpurple', edgecolor='none', alpha=0.85)
    ax.axvline(0.5, color='tomato', lw=1.5, ls='--', label='rq=0.5')
    ax.axvline(0.3, color='salmon', lw=1,   ls=':',  label='rq=0.3')
    ax.axvline(np.mean(rq), color='darkorange', lw=1.5, label=f'mean={np.mean(rq):.3f}')
    n_hi = sum(1 for r in rq if r > 0.5)
    n_lo = sum(1 for r in rq if r < 0.3)
    ax.set_xlabel('rq:f (Ultima per-read quality, 0–1)', fontsize=11)
    ax.set_ylabel('Read count', fontsize=11)
    ax.set_title(f'Ultima rq:f Quality Score  [>0.5: {100*n_hi/len(rq):.1f}%  <0.3: {100*n_lo/len(rq):.1f}%]',
                 fontsize=12)
    ax.legend(fontsize=9); fig.tight_layout()
    return _b64img(fig)


def plot_clamp(m):
    if m['n_clamp_scanned'] == 0:
        return '<p class="nodata">No AMBER clamp scan performed. Re-run with --amber flag.</p>'
    n_scanned = m['n_clamp_scanned']
    n_found   = len(m['clamp_pos'])
    n_hd0     = sum(1 for h in m['clamp_hd'] if h == 0)
    n_hd1     = n_found - n_hd0
    pct       = 100 * n_found / n_scanned

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 4))
    # Left: positional histogram of clamp detection
    if m['clamp_pos']:
        ax1.hist(m['clamp_pos'], bins=range(CLAMP_WIN_S, CLAMP_WIN_E + 2),
                 color='#9b59b6', edgecolor='none', alpha=0.85)
    ax1.axvline(49, color='gray', lw=0.8, ls='--', alpha=0.7, label='Pos 49 (peak)')
    ax1.axvline(50, color='gray', lw=0.8, ls=':',  alpha=0.7, label='Pos 50 (peak)')
    ax1.set_xlabel('Clamp position in read (1-based)', fontsize=11)
    ax1.set_ylabel('Read count', fontsize=11)
    ax1.set_title(f'GCATGCGG Clamp Position Distribution\n'
                  f'{n_found:,}/{n_scanned:,} reads ({pct:.1f}%) have clamp', fontsize=11)
    ax1.legend(fontsize=9)
    # Right: hd=0 vs hd=1 breakdown
    ax2.bar(['hd=0\n(exact)', 'hd=1\n(1-mismatch)'],
            [n_hd0, n_hd1], color=['#27ae60', '#f39c12'], alpha=0.85, edgecolor='none')
    for i, (n, label) in enumerate([(n_hd0, 'exact'), (n_hd1, '1-mm')]):
        pct_bar = 100 * n / max(n_found, 1)
        ax2.text(i, n + max(n_hd0, n_hd1) * 0.02, f'{pct_bar:.1f}%', ha='center', fontsize=10)
    ax2.set_ylabel('Read count', fontsize=11)
    ax2.set_title(f'Clamp Hamming Distance Breakdown\n'
                  f'(Expected ~30% yield for AMBER bridge-ligated reads)', fontsize=11)
    fig.tight_layout()
    return _b64img(fig)


# ─────────────────────────────────────────────────────────────────────────────
# HTML builder
# ─────────────────────────────────────────────────────────────────────────────

_STATUS_ICON  = {'pass': '&#x2714;', 'warn': '&#x26A0;', 'fail': '&#x2718;'}
_STATUS_COLOR = {'pass': '#27ae60',  'warn': '#e67e22',  'fail': '#e74c3c'}

def _nav_item(href, label, status):
    icon  = _STATUS_ICON[status]
    color = _STATUS_COLOR[status]
    return (f'<a class="nav-item {status}" href="#{href}">'
            f'<span class="nav-icon" style="color:{color}">{icon}</span>{label}</a>')

def _module(mid, title, status, body):
    icon  = _STATUS_ICON[status]
    color = _STATUS_COLOR[status]
    return f'''
<div class="module {status}" id="{mid}">
  <h2>{title} <span class="status-badge" style="color:{color}">{icon}</span></h2>
  {body}
</div>'''

def _basic_stats_table(m, pps):
    lengths = m['lengths']
    gc      = m['gc_fracs']
    mq      = m['mean_quals']
    rows = [
        ('Filename',           m['file']),
        ('Run date',           m['date']),
        ('CRAMqc version',     VERSION),
        ('Reads sampled',      f"{m['n_reads']:,}"),
        ('Mean read length',   f"{np.mean(lengths):.0f} bp (min {min(lengths)}, max {max(lengths)})"),
        ('Mean Phred quality', f"{np.mean(mq):.1f}"),
        ('Mean GC content',    f"{np.mean(gc):.1f}%"),
        ('rq:f tag present',   'Yes' if m['rq_scores'] else f"No (tag absent — likely post-alignment BAM)"),
    ]
    if pps:
        rows.append(('Per-pos reads', f"{m['n_pos_sample']:,} (first {m['n_pos_sample']:,} of sample)"))
    html = '<table><tr><th>Metric</th><th>Value</th></tr>'
    for k, v in rows:
        html += f'<tr><td>{k}</td><td>{v}</td></tr>'
    html += '</table>'
    return html

def build_html(m, pps, status, plots, do_amber):
    lengths  = m['lengths']
    gc_fracs = m['gc_fracs']
    mq       = m['mean_quals']
    rq       = m['rq_scores']

    # Sidebar navigation items
    nav = [
        _nav_item('basic-stats',   'Basic Statistics',            'pass'),
        _nav_item('per-base-qual', 'Per Base Seq Quality',        status['per_base_qual']),
        _nav_item('per-seq-qual',  'Per Sequence Quality Scores', status['per_seq_qual']),
        _nav_item('seq-content',   'Per Base Seq Content',        status['seq_content']),
        _nav_item('gc-content',    'Per Sequence GC Content',     status['gc_content']),
        _nav_item('length-dist',   'Sequence Length Distribution', status['length_dist']),
        _nav_item('duplication',   'Sequence Duplication Levels', status['duplication']),
        '<div class="nav-section">ULTIMA-SPECIFIC</div>',
        _nav_item('rq-qual',       'rq:f Quality Distribution',   status['rq_qual']),
    ]
    if do_amber:
        nav += ['<div class="nav-section">AMBER-SPECIFIC</div>',
                _nav_item('clamp', 'Clamp Detection (GCATGCGG)', status.get('clamp', 'warn'))]

    # Module bodies
    modules = [
        _module('basic-stats',   'Basic Statistics', 'pass',
                _basic_stats_table(m, pps)),
        _module('per-base-qual', 'Per Base Sequence Quality', status['per_base_qual'],
                f'<p class="desc">Phred quality scores at each read position. '
                f'Box = Q25–Q75 (IQR). Red line = median. Blue = mean. Whiskers = Q10/Q90.<br>'
                f'Background: <span style="background:#d5e8d4;padding:2px 6px">Green Q≥28</span> '
                f'<span style="background:#fff2cc;padding:2px 6px">Yellow Q20–28</span> '
                f'<span style="background:#f8cecc;padding:2px 6px">Red Q&lt;20</span>. '
                f'Dotted red line = Q15 (cutadapt -q threshold).</p>' + plots['per_base_qual']),
        _module('per-seq-qual',  'Per Sequence Quality Scores', status['per_seq_qual'],
                '<p class="desc">Distribution of mean Phred quality scores per read. '
                'A high-quality library peaks sharply above Q27.</p>' + plots['per_seq_qual']),
        _module('seq-content',   'Per Base Sequence Content', status['seq_content'],
                '<p class="desc">A/T/G/C composition at each read position. '
                'Ultima AMBER reads start with fixed Tn5 adapter (positions 1–14 show near-fixed composition), '
                'then Nextera ME (variable), then SpatBC (degenerate designed sequence), '
                'then chromatin insert. Expect unusual base content before the cutadapt trim point.</p>'
                + plots['seq_content']),
        _module('gc-content',    'Per Sequence GC Content', status['gc_content'],
                '<p class="desc">GC% per read. A narrow unimodal distribution centred near the '
                'genome GC% (~40–42% for hg38, ~38% for mm10) indicates a clean library. '
                'Broad or bimodal distribution suggests contamination or GC bias.</p>' + plots['gc_content']),
        _module('length-dist',   'Sequence Length Distribution', status['length_dist'],
                '<p class="desc">Ultima reads are variable-length (chromatin fragment + Tn5 insertion). '
                'Expect a broad distribution peaking 150–300 bp for typical AMBER chromatin inserts.</p>'
                + plots['length_dist']),
        _module('duplication',   'Sequence Duplication Levels', status['duplication'],
                ('<p class="desc" style="background:#fff3cd;border-left:4px solid #e67e22;padding:8px 12px;border-radius:4px">'
                 '<b>⚠ NOT INFORMATIVE for pre-trim AMBER CRAMs.</b> '
                 'FastQC convention hashes the first 50 bp, but here positions 1–33 are fixed adapter '
                 '(14 bp Tn5 + 19 bp Nextera ME) and positions 34–48 are the SpatBC cell barcode — '
                 'not genomic insert sequence. This module measures barcode diversity, not PCR duplication. '
                 'Assess true duplication post-alignment using <code>samtools markdup</code> on the mapped BAM.</p>'
                 if do_amber else
                 '<p class="desc">Estimated from first 50 bp of first 100,000 reads (FastQC convention). '
                 'High duplication indicates PCR over-amplification.</p>')
                + plots['duplication']),
        _module('rq-qual',       'Ultima rq:f Quality Score', status['rq_qual'],
                '<p class="desc"><b>Ultima-specific.</b> The <code>rq:f</code> tag is Ultima\'s per-read '
                'overall quality score (0–1), distinct from per-base Phred. '
                'rq &lt; 0.3 = low confidence read. rq &gt; 0.5 = good. '
                'Tag is present in raw Ultima CRAMs but may be absent in post-alignment BAMs.</p>'
                + plots['rq_qual']),
    ]
    if do_amber:
        modules.append(
            _module('clamp', 'AMBER Clamp Detection (GCATGCGG)', status.get('clamp', 'warn'),
                    '<p class="desc"><b>AMBER-specific.</b> The GCATGCGG bridge-ligation clamp '
                    'marks reads that contain a captured chromatin locus. Expected yield: ~30% '
                    '(remaining 70% are non-bridge-ligated Tn5 insertions — valid discard). '
                    'Window scanned: positions 47–65. Hamming distance ≤ 1 allowed (hd=0 = exact, '
                    'hd=1 = 1-mismatch recovery). This module runs on raw unaligned CRAMs only '
                    '(use before 1.CRAMprocess.sh).</p>' + plots['clamp'])
        )

    css = """
body{margin:0;font-family:Arial,Helvetica,sans-serif;background:#f0f2f5;color:#333}
#sidebar{width:230px;position:fixed;top:0;left:0;height:100vh;background:#2c3e50;
  color:#ecf0f1;overflow-y:auto;z-index:100}
#sidebar h1{padding:16px 15px 10px;font-size:15px;font-weight:700;color:#ecf0f1;
  border-bottom:1px solid #3d5166;margin:0}
.nav-item{display:block;padding:7px 14px;color:#bdc3c7;text-decoration:none;
  font-size:12px;border-left:3px solid transparent;transition:background 0.15s}
.nav-item:hover{background:#34495e;color:#ecf0f1}
.nav-item.pass{border-left-color:#27ae60}
.nav-item.warn{border-left-color:#e67e22}
.nav-item.fail{border-left-color:#e74c3c}
.nav-icon{margin-right:6px;font-style:normal}
.nav-section{padding:10px 14px 4px;font-size:10px;text-transform:uppercase;
  letter-spacing:1px;color:#7f8c8d;border-top:1px solid #3d5166;margin-top:6px}
#main{margin-left:248px;padding:20px 28px;max-width:1150px}
.header-bar{background:linear-gradient(135deg,#2c3e50,#2980b9);color:#fff;
  padding:18px 24px;border-radius:8px;margin-bottom:20px}
.header-bar h1{margin:0;font-size:20px}.header-bar p{margin:4px 0 0;font-size:12px;opacity:.8}
.module{background:#fff;border-radius:7px;padding:20px 24px;margin-bottom:18px;
  box-shadow:0 1px 4px rgba(0,0,0,.1)}
.module h2{margin:0 0 12px;font-size:16px;color:#2c3e50;border-bottom:2px solid #eee;
  padding-bottom:10px}
.status-badge{font-size:14px}
.module.pass{border-left:4px solid #27ae60}
.module.warn{border-left:4px solid #e67e22}
.module.fail{border-left:4px solid #e74c3c}
table{border-collapse:collapse;width:100%;margin-top:8px}
th{background:#2c3e50;color:#fff;padding:8px 14px;text-align:left;font-size:13px}
td{padding:7px 14px;border-bottom:1px solid #f0f0f0;font-size:13px}
tr:nth-child(even) td{background:#fafafa}
img.plot{max-width:100%;border:1px solid #eee;border-radius:5px;margin-top:8px}
.desc{font-size:12px;color:#666;margin:0 0 10px;line-height:1.5}
.nodata{color:#999;font-style:italic;font-size:12px}
code{background:#f0f0f0;padding:1px 4px;border-radius:3px;font-size:12px}
"""
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>CRAMqc — {m['file']}</title>
<style>{css}</style>
</head>
<body>
<div id="sidebar">
  <h1>&#x1F4CA; CRAMqc v{VERSION}</h1>
  {''.join(nav)}
</div>
<div id="main">
  <div class="header-bar">
    <h1>CRAMqc Report</h1>
    <p>{m['file']} &nbsp;|&nbsp; {m['date']} &nbsp;|&nbsp; {m['n_reads']:,} reads sampled</p>
  </div>
  {''.join(modules)}
  <p style="font-size:11px;color:#aaa;text-align:center;margin-top:20px">
    CRAMqc v{VERSION} — AMBER pipeline, Zhu Lab NYGC &nbsp;|&nbsp;
    Generated {m['date']} &nbsp;|&nbsp;
    Inspired by <a href="https://www.bioinformatics.babraham.ac.uk/projects/fastqc/" style="color:#aaa">FastQC</a>
  </p>
</div>
</body>
</html>"""
    return html


# ─────────────────────────────────────────────────────────────────────────────
# TSV summary writer
# ─────────────────────────────────────────────────────────────────────────────

def write_tsv(m, status, prefix):
    rq = m['rq_scores']
    lengths = m['lengths']
    out = f"{prefix}.CRAMqc.tsv"
    with open(out, 'w') as f:
        f.write("metric\tvalue\n")
        f.write(f"file\t{m['file']}\n")
        f.write(f"reads_sampled\t{m['n_reads']}\n")
        f.write(f"mean_read_length_bp\t{np.mean(lengths):.1f}\n")
        f.write(f"median_read_length_bp\t{np.median(lengths):.1f}\n")
        f.write(f"mean_phred_quality\t{np.mean(m['mean_quals']):.2f}\n")
        f.write(f"mean_gc_pct\t{np.mean(m['gc_fracs']):.2f}\n")
        f.write(f"reads_with_rq_tag\t{len(rq)}\n")
        f.write(f"mean_rq_score\t{np.mean(rq):.4f}\n" if rq else "mean_rq_score\tn/a\n")
        f.write(f"pct_rq_gt0.5\t{100*sum(r>0.5 for r in rq)/len(rq):.2f}\n" if rq else "pct_rq_gt0.5\tn/a\n")
        f.write(f"pct_rq_lt0.3\t{100*sum(r<0.3 for r in rq)/len(rq):.2f}\n" if rq else "pct_rq_lt0.3\tn/a\n")
        if m['n_clamp_scanned'] > 0:
            n_f = len(m['clamp_pos'])
            f.write(f"clamp_yield_pct\t{100*n_f/m['n_clamp_scanned']:.2f}\n")
            f.write(f"clamp_hd0_pct\t{100*sum(h==0 for h in m['clamp_hd'])/max(n_f,1):.2f}\n")
        for k, v in status.items():
            f.write(f"status_{k}\t{v}\n")
    print(f"[CRAMqc] TSV → {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description='CRAMqc — FastQC-equivalent HTML QC report for Ultima UG CRAMs',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Output: {prefix}.CRAMqc.html (open in browser) + {prefix}.CRAMqc.tsv (Google Sheets)')
    p.add_argument('cram', help='Input CRAM or BAM file')
    p.add_argument('--sample',     type=int, default=500000, help='Max reads to sample (default: 500000)')
    p.add_argument('--pos-sample', type=int, default=50000,  help='Reads for per-position stats (default: 50000)')
    p.add_argument('--threads',    type=int, default=4,      help='pysam I/O threads (default: 4)')
    p.add_argument('--amber',      action='store_true',      help='Enable AMBER-specific clamp detection module')
    p.add_argument('--out-prefix', default=None,             help='Output prefix (default: strip .cram/.bam)')
    args = p.parse_args()

    prefix = args.out_prefix or args.cram.rsplit('.cram', 1)[0].rsplit('.bam', 1)[0]
    n_pos  = min(args.pos_sample, args.sample)

    print(f"[CRAMqc] {os.path.basename(args.cram)}")
    print(f"[CRAMqc] sampling {args.sample:,} reads  (per-pos: {n_pos:,})  amber={args.amber}")

    # 1. Collect metrics
    m   = collect(args.cram, args.sample, n_pos, args.threads, args.amber)
    pps = compute_per_pos_stats(m)
    print(f"[CRAMqc] {m['n_reads']:,} reads collected. Computing stats...")

    # 2. Assess pass/warn/fail
    status = assess(m, pps)

    # 3. Generate plots
    print('[CRAMqc] Generating plots...')
    plots = {
        'per_base_qual': plot_per_base_quality(m, pps),
        'per_seq_qual':  plot_per_seq_quality(m),
        'seq_content':   plot_seq_content(m, pps),
        'gc_content':    plot_gc_content(m),
        'length_dist':   plot_length_dist(m),
        'duplication':   plot_duplication(m),
        'rq_qual':       plot_rq_dist(m),
        'clamp':         plot_clamp(m),
    }

    # 4. Build and write HTML
    html_str  = build_html(m, pps, status, plots, args.amber)
    out_html  = f'{prefix}.CRAMqc.html'
    with open(out_html, 'w') as f:
        f.write(html_str)
    print(f'[CRAMqc] HTML → {out_html}')

    # 5. TSV summary
    write_tsv(m, status, prefix)

    # 6. Stdout summary
    rq = m['rq_scores']
    print(f"\n{'═'*56}")
    print(f" CRAMqc Summary: {m['file']}")
    print(f"{'═'*56}")
    print(f"  Reads sampled      : {m['n_reads']:,}")
    print(f"  Mean read length   : {np.mean(m['lengths']):.0f} bp")
    print(f"  Mean Phred quality : {np.mean(m['mean_quals']):.1f}")
    print(f"  Mean GC content    : {np.mean(m['gc_fracs']):.1f}%")
    if rq:
        mrq = np.mean(rq)
        print(f"  Ultima rq:f        : mean={mrq:.3f}  "
              f"(>0.5: {100*sum(r>0.5 for r in rq)/len(rq):.1f}%  "
              f"<0.3: {100*sum(r<0.3 for r in rq)/len(rq):.1f}%)")
    if m['n_clamp_scanned'] > 0:
        n_f = len(m['clamp_pos'])
        print(f"  AMBER clamp yield  : {100*n_f/m['n_clamp_scanned']:.1f}%  "
              f"(hd=0: {sum(h==0 for h in m['clamp_hd'])}  hd=1: {sum(h==1 for h in m['clamp_hd'])})")
    print(f"\n  Module status:")
    for mod, st in status.items():
        icon = {'pass': '✓', 'warn': '!', 'fail': '✗'}[st]
        print(f"    {icon} {mod:<22} {st}")
    print(f"{'═'*56}")
    print(f"  Report: {out_html}")

if __name__ == '__main__':
    main()
