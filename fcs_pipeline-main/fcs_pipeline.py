"""
Demirer Lab — N. benthamiana Protoplast Flow Cytometry Pipeline
===============================================================
Sequential gating for transformation reporter / viability and target fluorescence detection.

Gate hierarchy
--------------
    P1         Intact protoplast scatter gate             FSC-A  vs  SSC-A
    P2         Singlet gate                               FSC-A  vs  FSC-Width
    P_MGOLD    Transformation reporter / viability gate   FSC-A  vs  525Blue-H
    MC         Target fluorescence readout gate           FSC-A  vs  610Yellow-H

Statistics
----------
For samples WITHOUT a transformation reporter signal (controls):
    Target% = MC ∩ P2  /  P2          (denominator = P2 singlets)

For samples WITH a transformation reporter signal (experimental constructs):
    Target% = MC ∩ P_MGOLD  /  P_MGOLD  (denominator = reporter⁺ cells)

Outputs
-------
    gating_strategy.svg              — 4-sample × 4-gate methods figure (SVG, dpi=600)
    final_fluorescence_gate.pdf   — channels-as-rows × samples-as-columns figure

Instrument  CytFLEX (Beckman Coulter)
Channels    525Blue-H  excitation 488 nm, emission ~530 nm  (FDA / mGold2t)
            610Yellow-H  excitation 561 nm, emission ~610 nm  (mCherry)

Dependencies
------------
    pip install fcsparser numpy matplotlib scipy

Author      Demirer Lab, Caltech
"""

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Imports
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
import fcsparser
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from scipy.ndimage import gaussian_filter
import warnings, os, copy

warnings.filterwarnings('ignore')


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# USER CONFIGURATION — edit this section for each experiment
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Directory containing .fcs files
DATA_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    '..', 'FCSfiles'
)

# Output directory (defaults to the same folder as this script)
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

# Channel names as stored in the FCS files
CH_525    = 'FDA 525Blue-H'       # transformation reporter / viability channel
CH_CHERRY = 'MCherry 610Yellow-H' # target fluorescence channel
FSC       = 'FSC-A'
SSC       = 'SSC-A'
FSC_W     = 'FSC-Width'

# Display names for the two fluorescence channels
# Change these to match your fluorophores (e.g. 'BFP', 'YFP')
REPORTER_LABEL = 'mGold2t'  # transformation reporter / viability marker (CH_525)
TARGET_LABEL   = 'mCherry'  # target fluorescence readout (CH_CHERRY)

# FSC axis scale factor: raw values are ~1e7, display as ×10⁵
S = 1 / 1e5

# Dot size multiplier: 1.0 = auto-scaled, 0.5 = half size (for dense samples)
DOT_SCALE = 0.5

# ── Samples for GATING STRATEGY figure ────────────────────────────────────────
# (filename, display_label, label_color)
# Rows are rendered top to bottom in this order.
GATING_SAMPLES = [
    ('wt Nb-noFDA.fcs',   'WT Nb\n(no FDA)',                    '#555555'),
    ('wt Nb-FDA.fcs',     'WT Nb\n(+FDA)',                      '#e07b22'),
    ('240 Nb no FDA.fcs', 'mCherry positive\ncontrol (no FDA)', '#1a7a6e'),
    ('240 Nb-FDA.fcs',    'mCherry positive\ncontrol (+FDA)',   '#2563ae'),
]

# ── Samples for SUPPLEMENTARY mCherry figure ──────────────────────────────────
# (filename, column_title, ch_525_name, is_control)
# is_control=True  → mCherry% denominator is P2 singlets
# is_control=False → mCherry% denominator is mGold2t⁺ cells
SUPP_SAMPLES = [
    ('wt Nb-noFDA.fcs',        'WT Nb',
     'FDA 525Blue-H',  True),
    ('240 Nb no FDA.fcs',      'mCherry positive\ncontrol',
     'FDA 525Blue-H',  True),
    ('536 mGold_low C gain.fcs',  'R2Tg intronized\n(mGold2t reporter)',
     'mGold 525Blue-H', False),
    ('502+536_low C gain2.fcs',   'R2Tg intronized (mGold2t)\n+ mCherry payload',
     'mGold 525Blue-H', False),
]

SUPP_SAMPLE_COLORS = ['#555555', '#1a7a6e', '#b8860b', '#8e44ad']


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# GATE COORDINATES
# All coordinates are in raw instrument units (FSC-A ~1–10 × 10⁶, etc.).
# Polygons are closed automatically in poly_gate() / draw_poly().
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

GATES = {
    # P1 — selects intact protoplast cloud on FSC-A vs SSC-A.
    # Excludes high-SSC aggregates and out-of-range debris.
    'P1': [
        (1_102_403,  1_200_000),
        (10_000_000, 3_800_000),
        (10_000_000, 1_700_000),
        (1_208_215,   130_199),
    ],

    # P2 — singlet gate on FSC-A vs FSC-Width.
    # Doublets appear above the diagonal as broad-width events.
    'P2': [
        (1_200_000,  2_500),
        (10_000_000, 8_600),
        (10_000_000, 6_500),
        (1_200_000,  1_200),
    ],

    # P_MGOLD — mGold2t⁺ / FDA⁺ gate on FSC-A vs 525Blue-H.
    # Lower FSC-A boundary set to 1 × 10⁶ to exclude small debris.
    # Upper boundary traces the autofluorescence shoulder to capture
    # live (FDA⁺) or mGold2t-expressing cells.
    'P_MGOLD': [
        (1_000_000,  10_000_000),
        (10_000_000, 10_000_000),
        (10_000_000,     22_000),
        (8_000_000,      22_000),
        (6_000_000,      22_000),
        (4_000_000,      21_000),
        (2_500_000,      18_000),
        (1_000_000,      14_000),
    ],

    # MC — mCherry⁺ gate on FSC-A vs 610Yellow-H.
    # Gate boundary set above chloroplast autofluorescence using the
    # mCherry positive control sample (240 Nb no FDA).
    'MC': [
        (1_000_000,  10_000_000),
        (10_000_000, 10_000_000),
        (10_000_000,      4_000),
        (1_000_000,       4_000),
    ],
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# STYLE
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# Pseudocolor scale: dark navy (sparse) → blue → cyan → green → yellow → red (dense)
_FLOW_COLORS = [
    (0.00, '#00003F'),
    (0.15, '#000080'),
    (0.35, '#0000FF'),
    (0.50, '#00CCFF'),
    (0.62, '#00FF80'),
    (0.72, '#FFFF00'),
    (0.84, '#FF6000'),
    (1.00, '#FF0000'),
]
CMAP_FLOW = LinearSegmentedColormap.from_list('flow_dark', _FLOW_COLORS, N=512)

GATE_KW = dict(color='black', lw=1.2, ls='-', zorder=5)

plt.rcParams.update({
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size':         9,
    'axes.labelsize':    9,
    'xtick.labelsize':   8,
    'ytick.labelsize':   8,
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'axes.linewidth':    0.7,
    'xtick.major.width': 0.7,
    'ytick.major.width': 0.7,
    'xtick.major.size':  3,
    'ytick.major.size':  3,
    'pdf.fonttype':      42,   # editable text in Illustrator
    'svg.fonttype':      'none',
})


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# CORE FUNCTIONS
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def poly_gate(d, cx, cy, poly):
    """Return boolean mask of events inside a polygon gate.

    Parameters
    ----------
    d    : pandas DataFrame — FCS data from fcsparser
    cx   : str — column name for x-axis channel
    cy   : str — column name for y-axis channel
    poly : list of (x, y) tuples — polygon vertices in raw instrument units

    Returns
    -------
    np.ndarray of bool, shape (n_events,)
    """
    pts = np.column_stack([d[cx].values, d[cy].values])
    return Path(poly + [poly[0]]).contains_points(pts)


def draw_poly(ax, poly, sx=1.0, sy=1.0, **kw):
    """Draw a closed polygon gate on an axes.

    Parameters
    ----------
    ax        : matplotlib Axes
    poly      : list of (x, y) tuples — polygon vertices
    sx, sy    : float — scale factors applied to x and y coordinates
                (use S = 1/1e5 for the FSC axis when displaying as ×10⁵)
    **kw      : passed to ax.plot (e.g. color, lw, ls)
    """
    xs = [v[0] * sx for v in poly] + [poly[0][0] * sx]
    ys = [v[1] * sy for v in poly] + [poly[0][1] * sy]
    ax.plot(xs, ys, **kw)


def pcolor_density(ax, x, y, log_y=False, s=None, subsample=None):
    """Pseudocolor density scatter plot.

    Events are colored by local density using a 2-D histogram smoothed
    with a Gaussian kernel, then sorted so dense events are plotted last
    (on top). The colormap runs dark navy (sparse) → red (dense).

    Parameters
    ----------
    ax        : matplotlib Axes
    x, y      : array-like — coordinates to plot
    log_y     : bool — if True, y values are log10-transformed before
                binning and the caller should set ax.set_yscale('log')
    s         : float or None — marker size; if None, auto-scaled by
                sqrt(n_events) × DOT_SCALE
    subsample : int or None — if given, randomly subsample to this many
                events before plotting (useful for very large datasets;
                None = plot all events)
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) == 0:
        return
    if log_y:
        y = np.clip(y, 100, None)
    if subsample is not None and len(x) > subsample:
        rng = np.random.default_rng(42)
        idx = rng.choice(len(x), subsample, replace=False)
        x, y = x[idx], y[idx]
    if s is None:
        s = float(np.clip(800 / np.sqrt(max(len(x), 1)), 2.0, 8.0)) * DOT_SCALE
    y_h   = np.log10(y) if log_y else y
    bins  = min(300, max(80, len(x) // 15))
    H, xe, ye = np.histogram2d(x, y_h, bins=bins)
    H_s   = gaussian_filter(H.astype(float), sigma=1.0)
    ix    = np.clip(np.searchsorted(xe[1:-1], x),   0, H.shape[0] - 1)
    iy    = np.clip(np.searchsorted(ye[1:-1], y_h), 0, H.shape[1] - 1)
    z     = H_s[ix, iy]
    order = np.argsort(z)
    ax.scatter(x[order], y[order],
               c=z[order] + 1,
               cmap=CMAP_FLOW, norm=LogNorm(),
               s=s, linewidths=0, rasterized=True)


def style_log_axis(ax, axis='y'):
    """Apply log-scale tick formatting (major powers of 10, minor ticks)."""
    loc    = mticker.LogLocator(base=10, numticks=6)
    minloc = mticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=50)
    fmt    = mticker.LogFormatterMathtext(base=10, labelOnlyBase=True)
    if axis == 'y':
        ax.yaxis.set_major_locator(loc)
        ax.yaxis.set_major_formatter(fmt)
        ax.yaxis.set_minor_locator(minloc)
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        ax.tick_params(axis='y', which='minor', length=2, width=0.4)


def stat_label(ax, num, denom, fontsize=8):
    """Print 'n=num/denom\\npct%' in the bottom-right corner of an axes."""
    pct = 100 * num / denom if denom else 0
    ax.text(0.97, 0.04, f'n={num:,}/{denom:,}\n{pct:.1f}%',
            fontsize=fontsize, fontweight='bold', color='black',
            ha='right', va='bottom', transform=ax.transAxes)


def load_and_gate(fname, ch_525=None, data_dir=None):
    """Load one FCS file and apply the full sequential gate hierarchy.

    Parameters
    ----------
    fname    : str — filename (relative to data_dir)
    ch_525   : str or None — 525Blue-H channel name; defaults to CH_525
    data_dir : str or None — directory containing FCS files; defaults to DATA_DIR

    Returns
    -------
    dict with keys:
        d            pandas DataFrame (all events)
        g            gate dict (deep copy of GATES)
        p1, p2       boolean masks for P1 and P2 gates
        p_mg         P_MGOLD ∩ P2 mask (mGold2t⁺ / FDA⁺)
        p_mc         MC ∩ P_MGOLD mask (mCherry⁺ within mGold2t⁺)
        p_mc_from_p2 MC ∩ P2 mask (mCherry⁺ relative to P2 — for controls)
        N, n1, n2, n_mg, n_mc, n_mc_from_p2   event counts
    """
    if ch_525 is None:
        ch_525 = CH_525
    if data_dir is None:
        data_dir = DATA_DIR
    _, d    = fcsparser.parse(os.path.join(data_dir, fname), reformat_meta=True)
    g       = copy.deepcopy(GATES)
    p1      = poly_gate(d, FSC, SSC,   g['P1'])
    p2      = poly_gate(d, FSC, FSC_W, g['P2'])       & p1
    p_mg    = poly_gate(d, FSC, ch_525,     g['P_MGOLD']) & p2
    p_mc    = poly_gate(d, FSC, CH_CHERRY,  g['MC'])      & p_mg
    p_mc_p2 = poly_gate(d, FSC, CH_CHERRY,  g['MC'])      & p2
    return dict(
        d=d, g=g,
        p1=p1, p2=p2, p_mg=p_mg, p_mc=p_mc, p_mc_from_p2=p_mc_p2,
        N=len(d), n1=p1.sum(), n2=p2.sum(),
        n_mg=p_mg.sum(), n_mc=p_mc.sum(), n_mc_from_p2=p_mc_p2.sum(),
    )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# FIGURE 1 — GATING STRATEGY
# 4 samples (rows) × 4 gate steps (columns), with annotation text below.
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def make_gating_strategy_figure(samples=None, out_dir=None):
    """Generate and save the gating strategy figure.

    Parameters
    ----------
    samples : list of (fname, label, color) — defaults to GATING_SAMPLES
    out_dir : str — output directory; defaults to OUT_DIR
    """
    if samples is None:
        samples = GATING_SAMPLES
    if out_dir is None:
        out_dir = OUT_DIR

    print("Gating strategy figure — loading samples (no subsampling)…")
    rows = []
    for fname, label, color in samples:
        r = load_and_gate(fname)
        r['label'] = label
        r['color'] = color
        rows.append(r)
        print(f"  {label.replace(chr(10), ' ')}: "
              f"N={r['N']:,}  P1={r['n1']:,}  P2={r['n2']:,}  "
              f"{REPORTER_LABEL}⁺={r['n_mg']:,}  {TARGET_LABEL}⁺={r['n_mc']:,}")

    # ── Shared axis limits ────────────────────────────────────────────────────
    xlim     = (-2, 120)
    ylim_ssc = (-2, 60)

    fw_vals     = np.concatenate([r['d'].loc[r['p1'], FSC_W].values for r in rows])
    ylim_fw     = (0, max(np.percentile(fw_vals, 99.99),
                          max(v[1] for v in GATES['P2']) * 1.5))

    fda_vals    = np.concatenate([np.clip(r['d'].loc[r['p2'], CH_525].values, 100, None)
                                  for r in rows])
    ylim_525    = (
        10 ** (np.floor(np.log10(np.percentile(fda_vals,  0.1))) - 0.15),
        max(10 ** (np.ceil(np.log10(np.percentile(fda_vals, 99.9))) + 0.3),
            max(v[1] for v in GATES['P_MGOLD']) * 3.0)
    )

    mc_vals     = np.concatenate([np.clip(r['d'].loc[r['p2'], CH_CHERRY].values, 100, None)
                                  for r in rows])
    ylim_mc     = (
        10 ** (np.floor(np.log10(np.percentile(mc_vals,  0.1))) - 0.15),
        max(10 ** (np.ceil(np.log10(np.percentile(mc_vals, 99.9))) + 0.3),
            max(v[1] for v in GATES['MC']) * 3.0)
    )

    # ── Layout ────────────────────────────────────────────────────────────────
    ndata  = len(rows)
    ncols  = 4
    fig    = plt.figure(figsize=(15, 2.6 * ndata + 3.8), dpi=300)
    gs_fig = gridspec.GridSpec(
        ndata + 1, ncols,
        height_ratios=[2.6] * ndata + [2.8],
        hspace=0.06, wspace=0.30,
        left=0.09, right=0.98, top=0.93, bottom=0.03,
    )
    axes    = np.array([[fig.add_subplot(gs_fig[ri, ci])
                         for ci in range(ncols)]
                        for ri in range(ndata)])
    ann_axs = [fig.add_subplot(gs_fig[ndata, ci]) for ci in range(ncols)]
    for ax in ann_axs:
        ax.axis('off')

    # Column headers
    col_headers = ['All events', 'Gate: cells (P1)',
                   'Gate: singlets (P2)', f'Gate: {REPORTER_LABEL}⁺']
    for ci, hdr in enumerate(col_headers):
        axes[0, ci].set_title(hdr, fontsize=10, fontweight='bold', pad=8)

    # ── Data panels ───────────────────────────────────────────────────────────
    for ri, r in enumerate(rows):
        d, g  = r['d'], r['g']
        p1, p2, p_mg, p_mc = r['p1'], r['p2'], r['p_mg'], r['p_mc']
        N, n1, n2, n_mg, n_mc = r['N'], r['n1'], r['n2'], r['n_mg'], r['n_mc']
        is_bot = (ri == ndata - 1)

        # Col 0 — All events: FSC-A vs SSC-A
        ax = axes[ri, 0]
        pcolor_density(ax, d[FSC].values * S, d[SSC].values * S)
        draw_poly(ax, g['P1'], S, S, **GATE_KW)
        stat_label(ax, n1, N)
        ax.set(xlim=xlim, ylim=ylim_ssc)
        ax.set_ylabel('SSC-A (×10⁵)')
        if is_bot: ax.set_xlabel('FSC-A (×10⁵)')
        else:      ax.tick_params(labelbottom=False)

        # Col 1 — P1 events: FSC-A vs FSC-Width
        ax = axes[ri, 1]
        pcolor_density(ax, d.loc[p1, FSC].values * S, d.loc[p1, FSC_W].values)
        draw_poly(ax, g['P2'], S, 1.0, **GATE_KW)
        stat_label(ax, n2, n1)
        ax.set(xlim=xlim, ylim=ylim_fw)
        ax.set_ylabel('FSC-Width')
        if is_bot: ax.set_xlabel('FSC-A (×10⁵)')
        else:      ax.tick_params(labelbottom=False)

        # Col 2 — P2 singlets: FSC-A vs 525Blue-H (FDA / mGold2t)
        ax = axes[ri, 2]
        pcolor_density(ax, d.loc[p2, FSC].values * S,
                       np.clip(d.loc[p2, CH_525].values, 100, None),
                       log_y=True)
        draw_poly(ax, g['P_MGOLD'], S, 1.0, **GATE_KW)
        stat_label(ax, n_mg, n2)
        ax.set_yscale('log')
        ax.set(xlim=xlim, ylim=ylim_525)
        style_log_axis(ax)
        ax.set_ylabel(f'{REPORTER_LABEL}  525Blue-H')
        if is_bot: ax.set_xlabel('FSC-A (×10⁵)')
        else:      ax.tick_params(labelbottom=False)

        # Col 3 — FDA⁺/mGold2t⁺ events: FSC-A vs mCherry 610Yellow-H
        ax = axes[ri, 3]
        if n_mg > 0:
            pcolor_density(ax, d.loc[p_mg, FSC].values * S,
                           np.clip(d.loc[p_mg, CH_CHERRY].values, 100, None),
                           log_y=True)
        draw_poly(ax, g['MC'], S, 1.0, **GATE_KW)
        stat_label(ax, n_mc, n_mg)
        ax.set_yscale('log')
        ax.set(xlim=xlim, ylim=ylim_mc)
        style_log_axis(ax)
        ax.set_ylabel(f'{TARGET_LABEL}  610Yellow-H')
        if is_bot: ax.set_xlabel('FSC-A (×10⁵)')
        else:      ax.tick_params(labelbottom=False)

    # ── Annotation text boxes ─────────────────────────────────────────────────
    col_annotations = [
        ("All acquired events shown as FSC-A (cell size) vs SSC-A (granularity).\n"
         "Gate P1 selects the main intact protoplast population, excluding\n"
         "high-SSC aggregates. High SSC-A events may represent high-granularity\n"
         "mesophyll cells and can be analysed separately if needed."),
        ("P1-gated events shown as FSC-A vs FSC-Width (pulse width). Singlets\n"
         "form a tight diagonal band. Gate P2 excludes doublets and multiplets,\n"
         "which appear above the diagonal as broad-width events. The majority\n"
         "of P1-gated events fall within P2 (minimal doublet formation)."),
        (f"P2 singlets shown as FSC-A vs {REPORTER_LABEL} 525Blue-H. Gate selects\n"
         f"transformation reporter-positive or viable cells ({REPORTER_LABEL}⁺).\n"
         f"Adjust the P_MGOLD gate coordinates if your reporter population\n"
         f"shifts on your instrument or under different experimental conditions."),
        (f"{REPORTER_LABEL}⁺ singlets shown as FSC-A vs {TARGET_LABEL} 610Yellow-H.\n"
         f"Gate boundary is set above autofluorescence using a {TARGET_LABEL}\n"
         f"positive control sample. Event counts and percentages within each\n"
         f"gate are shown in the lower right corner of each panel."),
    ]
    box_style = dict(boxstyle='round,pad=0.5', facecolor='#f4f4f4',
                     edgecolor='#bbbbbb', linewidth=0.8)
    for ax, txt in zip(ann_axs, col_annotations):
        ax.text(0.5, 0.88, txt, ha='center', va='top',
                transform=ax.transAxes, fontsize=7.5,
                color='#222222', linespacing=1.45, bbox=box_style)

    # ── Arrows between column headers ─────────────────────────────────────────
    fig.canvas.draw()
    for ci in range(ncols - 1):
        pos_l = axes[0, ci].get_position()
        pos_r = axes[0, ci + 1].get_position()
        y_arr = pos_l.y1 + 0.030
        fig.add_artist(mpatches.FancyArrowPatch(
            posA=(pos_l.x1 + 0.004, y_arr),
            posB=(pos_r.x0 - 0.004, y_arr),
            transform=fig.transFigure,
            arrowstyle='->', color='#666666', lw=1.3,
            mutation_scale=12, zorder=10,
        ))

    # ── Sample labels (left margin) ───────────────────────────────────────────
    for ri, r in enumerate(rows):
        pos = axes[ri, 0].get_position()
        y_c = (pos.y0 + pos.y1) / 2
        fig.text(0.010, y_c, r['label'],
                 ha='center', va='center', rotation=90,
                 fontsize=9, fontweight='bold', color=r['color'],
                 linespacing=1.4)

    # ── Save ──────────────────────────────────────────────────────────────────
    path = os.path.join(out_dir, 'gating_strategy.svg')
    fig.savefig(path, format='svg', dpi=600, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# FIGURE 2 — FINAL FLUORESCENCE GATE FIGURE
# Channels as rows (reporter top, target bottom) × samples as columns.
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

def make_supp_mcherry_figure(samples=None, out_dir=None, filename='final_fluorescence_gate.pdf'):
    """Generate and save the supplementary mGold2t / mCherry figure.

    Parameters
    ----------
    samples  : list of (fname, col_title, ch_525, is_control)
               defaults to SUPP_SAMPLES
    out_dir  : str — output directory; defaults to OUT_DIR
    filename : str — output filename (PDF)
    """
    if samples is None:
        samples = SUPP_SAMPLES
    if out_dir is None:
        out_dir = OUT_DIR

    print("Final fluorescence gate figure — loading samples…")
    cols = []
    for fname, col_title, ch_525, is_control in samples:
        r = load_and_gate(fname, ch_525=ch_525)
        r['title']      = col_title
        r['ch_525']     = ch_525
        r['is_control'] = is_control
        cols.append(r)
        print(f"  {col_title.replace(chr(10), ' ')}: "
              f"N={r['N']:,}  P2={r['n2']:,}  "
              f"{REPORTER_LABEL}⁺={r['n_mg']:,}  "
              f"{TARGET_LABEL}⁺={r['n_mc']:,}  {TARGET_LABEL}⁺(from P2)={r['n_mc_from_p2']:,}")

    # ── Shared axis limits ────────────────────────────────────────────────────
    xlim = (-2, 120)

    b525_all = np.concatenate([
        np.clip(c['d'].loc[c['p2'], c['ch_525']].values, 100, None)
        for c in cols])
    ylim_525 = (
        10 ** (np.floor(np.log10(np.percentile(b525_all,  0.1))) - 0.15),
        max(10 ** (np.ceil(np.log10(np.percentile(b525_all, 99.9))) + 0.3),
            max(v[1] for v in GATES['P_MGOLD']) * 3.0)
    )

    mc_masks  = [c['p2'] if c['is_control'] else c['p_mg'] for c in cols]
    mc_all    = np.concatenate([
        np.clip(c['d'].loc[m, CH_CHERRY].values, 100, None)
        for c, m in zip(cols, mc_masks) if m.sum() > 0])
    ylim_mc   = (
        10 ** (np.floor(np.log10(np.percentile(mc_all,  0.1))) - 0.15),
        max(10 ** (np.ceil(np.log10(np.percentile(mc_all, 99.9))) + 0.3),
            max(v[1] for v in GATES['MC']) * 3.0)
    )

    # ── Layout ────────────────────────────────────────────────────────────────
    ncols_fig = len(cols)
    fig, axes = plt.subplots(2, ncols_fig,
                             figsize=(4.0 * ncols_fig, 7.6),
                             dpi=300)
    fig.subplots_adjust(left=0.10, right=0.99, top=0.88,
                        bottom=0.09, wspace=0.28, hspace=0.18)

    row_labels = [f'{REPORTER_LABEL}\n(525Blue-H)', f'{TARGET_LABEL}\n(610Yellow-H)']
    row_colors = ['#b8860b', '#b03060']

    for ci, (c, scol) in enumerate(zip(cols, SUPP_SAMPLE_COLORS)):
        axes[0, ci].set_title(c['title'],
                              fontsize=9, fontweight='bold',
                              color=scol, pad=14, linespacing=1.4)

    for ci, c in enumerate(cols):
        d, g   = c['d'], c['g']
        p2     = c['p2']
        p_mg   = c['p_mg']

        # Row 0 — mGold2t / FDA (525Blue-H), always gated from P2
        ax = axes[0, ci]
        pcolor_density(ax,
                       d.loc[p2, FSC].values * S,
                       np.clip(d.loc[p2, c['ch_525']].values, 100, None),
                       log_y=True)
        draw_poly(ax, g['P_MGOLD'], S, 1.0, **GATE_KW)
        pct = 100 * c['n_mg'] / c['n2'] if c['n2'] else 0
        stat_label(ax, c['n_mg'], c['n2'], fontsize=9)
        ax.text(0.50, 1.01, 'P2 singlets',
                fontsize=8, color='#444444', style='italic',
                ha='center', va='bottom', transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set(xlim=xlim, ylim=ylim_525)
        style_log_axis(ax)
        ax.set_ylabel(f'{REPORTER_LABEL} 525Blue-H')
        ax.tick_params(axis='x', labelbottom=False)
        if ci > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)

        # Row 1 — mCherry (610Yellow-H)
        # Controls:   plot P2 singlets, denominator = n2
        # Constructs: plot mGold2t⁺ events, denominator = n_mg
        ax       = axes[1, ci]
        mc_mask  = p2 if c['is_control'] else p_mg
        pop_lbl  = 'P2 singlets' if c['is_control'] else f'{REPORTER_LABEL}⁺ cells'

        if c['is_control']:
            n_num   = c['n_mc_from_p2']
            n_denom = c['n2']
        else:
            n_num   = c['n_mc']
            n_denom = c['n_mg']

        if mc_mask.sum() > 0:
            pcolor_density(ax,
                           d.loc[mc_mask, FSC].values * S,
                           np.clip(d.loc[mc_mask, CH_CHERRY].values, 100, None),
                           log_y=True)
        draw_poly(ax, g['MC'], S, 1.0, **GATE_KW)
        stat_label(ax, n_num, n_denom, fontsize=9)
        ax.text(0.50, 1.01, pop_lbl,
                fontsize=8, color='#444444', style='italic',
                ha='center', va='bottom', transform=ax.transAxes)
        ax.set_yscale('log')
        ax.set(xlim=xlim, ylim=ylim_mc)
        style_log_axis(ax)
        ax.set_ylabel(f'{TARGET_LABEL} 610Yellow-H')
        ax.set_xlabel('FSC-A (×10⁵)')
        if ci > 0:
            ax.set_ylabel('')
            ax.tick_params(labelleft=False)

    # ── Row labels (left margin) ───────────────────────────────────────────────
    fig.canvas.draw()
    for ri, (rlabel, rcol) in enumerate(zip(row_labels, row_colors)):
        pos = axes[ri, 0].get_position()
        y_c = (pos.y0 + pos.y1) / 2
        fig.text(0.015, y_c, rlabel,
                 ha='center', va='center', rotation=90,
                 fontsize=9, fontweight='bold', color=rcol, linespacing=1.4)

    # ── Save ──────────────────────────────────────────────────────────────────
    path = os.path.join(out_dir, filename)
    fig.savefig(path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {path}")
    plt.close(fig)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# RUN
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if __name__ == '__main__':
    make_gating_strategy_figure()
    make_supp_mcherry_figure()
    print("\nAll figures complete.")
