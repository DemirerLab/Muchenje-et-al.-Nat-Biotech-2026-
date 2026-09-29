#!/usr/bin/env python3
"""
junction_seq_table.py
─────────────────────
Demirer Lab — R2 3' junction PCR nanopore read analysis.

Orients reads using PCR primers, extracts a fixed window around the
junction (3' UTR payload → 25S rDNA), groups identical junction sequences,
and renders a publication-quality color-coded nucleotide table ranked by
frequency.

No nick site classification is performed. The raw junction sequences are
shown as observed.

Usage
-----
    python junction_seq_table.py

Edit the Config section below to point at your FASTQ file and set output
paths and window sizes.

Dependencies
------------
    matplotlib >= 3.5
    numpy (for nothing critical; can be removed if desired)

Authors
-------
    Demirer Lab, Caltech CCE
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import Counter

# ─── Config ───────────────────────────────────────────────────────────────────

FASTQ  = '/Users/doram/Desktop/2026 spring/20260507_Athaliana_seq/Athaliana_3_junction_PCRresults/54QHSZ_fastq/54QHSZ_2_430+431.B.fastq'
OUTDIR = '/Users/doram/Desktop/2026 spring/Athaliana_3_junctionPCR_seq'
LABEL  = '54QHSZ_430+431B'

# Template sequence (full amplicon, used to define the nick position anchor)
TEMPLATE = (
    'AGTTAGGTAGCGGATAGTAGGTAGGAACAGACTTTTACTATTTCATAACGCGTCAATTACCACCTGATTTGGACC'
    'AATTCACGGGATTTGTCCAAGGTGGACGGGCCACCTTTACTTAACCCGGAAAAGGAACATATATAATTTATGTGTG'
    'TTCGATAAATAGCCAAATGCCTCGTCATCTAATTAGTGACGCGCATGAATGGATTAACGAGATTCCCACTGTCCCT'
    'GTCTACTATCCAGCGAAACCACAGCCAAGGGAACGGGCTTGGCAGAATCAGCGGGGAAAGAAGACCCTGTTGAGCTT'
    'GACTCTAGTCCGACTTTGTGAAATGACTTGAGAGGTGTAGGATAAGTGGGAGCTTCGGCGCAAGTGAAATACCACT'
    'ACTTTTAACGTTATTTTACTTACTCCGTGAATCGGAGGCGGGG'
)

# Nick position: index in TEMPLATE where the 25S rDNA sequence begins
NICK_POS = 160

# PCR primers
FWD_PRIMER = 'AGTTAGGTAGCGGATAGTAGGTAGGAA'
REV_PRIMER = 'CCCCGCCTCCGATTCACGG'

# Window to display on each side of the nick
UP_WIN   = 11   # bp from 3' UTR (payload) side
DOWN_WIN = 11   # bp into 25S rDNA

# How many top junction sequences to show in the table
TOP_N = 12

# Primer mismatch tolerance
MAX_MM_PRIMER = 2
MAX_MM_ANCHOR = 3

# Nucleotide colors
NT_COLORS = {
    'A': '#F0E890',   # pale yellow
    'T': '#B8A8D8',   # lavender
    'G': '#9CC89C',   # sage green
    'C': '#90B8E0',   # light blue
    '-': '#E0E0E0',   # gap / missing
    '?': '#F0F0F0',   # unknown
}

# ─── Utilities ────────────────────────────────────────────────────────────────

RC_TABLE = str.maketrans('ACGTacgt', 'TGCAtgca')

def rc(seq):
    return seq.translate(RC_TABLE)[::-1]


def find_primer(seq, primer, max_mm=MAX_MM_PRIMER):
    """Sliding-window primer match allowing up to max_mm mismatches."""
    plen = len(primer)
    best_pos, best_mm = -1, max_mm + 1
    for i in range(len(seq) - plen + 1):
        mm = sum(a != b for a, b in zip(primer, seq[i:i + plen]))
        if mm < best_mm:
            best_mm, best_pos = mm, i
    return best_pos if best_mm <= max_mm else -1


def orient(seq):
    """Return read in forward orientation (primer trimmed) or (None, '?')."""
    rev_rc = rc(REV_PRIMER)
    fp = find_primer(seq, FWD_PRIMER)
    if fp >= 0:
        return seq[fp:], '+'
    fp = find_primer(seq, rev_rc)
    if fp >= 0:
        return seq[fp:], '+'
    seq2 = rc(seq)
    fp = find_primer(seq2, FWD_PRIMER)
    if fp >= 0:
        return seq2[fp:], '+'
    return None, '?'


def find_anchor(seq, anchor, max_mm=MAX_MM_ANCHOR):
    """Sliding-window anchor match allowing up to max_mm mismatches."""
    plen = len(anchor)
    best_pos, best_mm = -1, max_mm + 1
    for i in range(len(seq) - plen + 1):
        mm = sum(a != b for a, b in zip(anchor, seq[i:i + plen]))
        if mm < best_mm:
            best_mm, best_pos = mm, i
    return best_pos if best_mm <= max_mm else -1

# ─── FASTQ parsing ────────────────────────────────────────────────────────────

def parse_fastq(path):
    reads = []
    with open(path) as f:
        while True:
            name = f.readline().strip()
            if not name:
                break
            seq  = f.readline().strip()
            f.readline()   # +
            f.readline()   # quality
            reads.append(seq)
    return reads

# ─── Junction extraction ──────────────────────────────────────────────────────

# 20 bp anchor ending exactly at the nick site
ANCHOR_UP = TEMPLATE[NICK_POS - 20 : NICK_POS]

def extract_junction(seq):
    """
    Locate the nick site in an oriented read and return the flanking window.

    Returns
    -------
    (upstream_11bp, downstream_11bp) : (str, str)
        Strings of exactly UP_WIN and DOWN_WIN characters.
        Positions missing from the read are filled with '-'.
    Returns None if the anchor is not found.
    """
    pos = find_anchor(seq, ANCHOR_UP)
    if pos < 0:
        return None

    nick = pos + len(ANCHOR_UP)   # index of first rDNA base in read

    # Upstream window (3' UTR side)
    up_start = nick - UP_WIN
    if up_start < 0:
        upstream = '-' * (-up_start) + seq[0:nick]
    else:
        upstream = seq[up_start:nick]

    # Downstream window (25S rDNA side)
    dn_end = nick + DOWN_WIN
    if dn_end > len(seq):
        downstream = seq[nick:] + '-' * (dn_end - len(seq))
    else:
        downstream = seq[nick:dn_end]

    return upstream, downstream

# ─── Main ─────────────────────────────────────────────────────────────────────

reads = parse_fastq(FASTQ)
total = len(reads)
print(f'Total reads   : {total}')

n_oriented  = 0
n_extracted = 0
junctions   = []

for seq in reads:
    oriented, strand = orient(seq)
    if strand == '?':
        continue
    n_oriented += 1

    result = extract_junction(oriented)
    if result is None:
        continue
    n_extracted += 1
    junctions.append(result)   # (upstream_str, downstream_str)

print(f'Oriented      : {n_oriented}')
print(f'Junction found: {n_extracted}')

# Count unique junction sequences
counter = Counter(junctions)
top     = counter.most_common(TOP_N)

print(f'\nTop {TOP_N} junction sequences (% of {n_extracted} extracted reads):')
for (up, dn), cnt in top:
    pct = cnt / n_extracted * 100
    print(f'  {up}|{dn}  {cnt} reads  {pct:.2f}%')

# ─── Figure ───────────────────────────────────────────────────────────────────

REF_UP   = TEMPLATE[NICK_POS - UP_WIN : NICK_POS]
REF_DOWN = TEMPLATE[NICK_POS : NICK_POS + DOWN_WIN]

# Build row list: reference first, then top patterns
rows = [{'up': REF_UP, 'dn': REF_DOWN, 'count': None, 'pct': None, 'is_ref': True}]
for (up, dn), cnt in top:
    rows.append({
        'up':    up,
        'dn':    dn,
        'count': cnt,
        'pct':   cnt / n_extracted * 100,
        'is_ref': False,
    })

n_rows   = len(rows)
n_cols   = UP_WIN + DOWN_WIN   # no gap column — just the raw window

CELL_W = 0.55
CELL_H = 0.50
fig_w  = n_cols * CELL_W + 7.0
fig_h  = n_rows * CELL_H + 2.5

plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams.update({'pdf.fonttype': 42, 'svg.fonttype': 'none'})

fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=200)
ax.set_xlim(0, n_cols)
ax.set_ylim(0, n_rows)
ax.axis('off')

RIGHT_X = n_cols + 0.3

for ri, row in enumerate(rows):
    yi   = n_rows - ri - 1   # row 0 at top
    seq  = list(row['up']) + list(row['dn'])

    # Right-side label
    if row['is_ref']:
        ax.text(RIGHT_X, yi + 0.5, 'Reference',
                ha='left', va='center', fontsize=18, fontweight='bold',
                color='#333333')
    else:
        ax.text(RIGHT_X, yi + 0.5,
                f"{row['pct']:.1f}%  ({row['count']} reads)",
                ha='left', va='center', fontsize=17, color='#555555')

    # Nucleotide cells
    for ci, nt in enumerate(seq):
        color = NT_COLORS.get(nt.upper(), NT_COLORS['?'])
        rect  = plt.Rectangle((ci, yi), 1, 1,
                               facecolor=color, edgecolor='white', linewidth=1.5)
        ax.add_patch(rect)
        ax.text(ci + 0.5, yi + 0.5, nt.upper(),
                ha='center', va='center',
                fontsize=16, color='#222222')

    # Blue border on reference row
    if row['is_ref']:
        border = plt.Rectangle((0, yi), n_cols, 1,
                                facecolor='none', edgecolor='#222288', linewidth=2.5)
        ax.add_patch(border)

# Dashed divider at nick site
ax.axvline(UP_WIN, color='#333333', linewidth=1.5, linestyle='--', zorder=5)

# Region labels below
ax.text(UP_WIN / 2,       -0.35, "3' UTR (payload)",
        ha='center', va='top', fontsize=17, fontweight='bold', color='#555555')
ax.text(UP_WIN + DOWN_WIN / 2, -0.35, '25S rDNA',
        ha='center', va='top', fontsize=17, fontweight='bold', color='#555555')

fig.suptitle(f"{LABEL} — R2 3' junction sequences",
             fontsize=20, fontweight='bold', y=1.01)

os.makedirs(OUTDIR, exist_ok=True)
for ext in ('png', 'svg'):
    p = os.path.join(OUTDIR, f'{LABEL}_junction_table.{ext}')
    fig.savefig(p, dpi=200, bbox_inches='tight')
    print(f'Saved: {p}')
plt.close()

print('\nDone.')
