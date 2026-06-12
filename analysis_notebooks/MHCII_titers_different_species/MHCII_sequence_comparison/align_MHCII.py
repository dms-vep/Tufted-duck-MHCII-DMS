"""Align MHCII alpha chain sequences and visualize the binding-site regions.

Everything is driven by config.yaml:

    MHCII_sequences                input FASTA of MHCII alpha chain sequences
    compared_to                    reference sequence (top row; residue numbering)
    Tufted_duck_MHCII_binding_sites  regions in reference numbering, e.g.
                                     "31-39, 44-47, 59-66, 81-96"
    exclude                        sequence name(s) to drop entirely
                                     (string, comma list, or YAML list)
    Tufted_duck_MHCII_binding_data   (optional) CSV of deep-mutational-scanning
                                     data with columns 'ORF_site', 'site',
                                     'mutant', and 'H5 HA binding escape'.
                                     'ORF_site' is the integer position in the
                                     reference numbering used for matching;
                                     'site' is the display label shown on the
                                     figure. When given, alignment squares are
                                     colored by the HA binding escape of the
                                     residue each sequence carries at that site
                                     (tufted-duck wildtype = 0).

Usage:
    python align_MHCII.py [config.yaml]

Outputs:
    MHCII_aligned.fasta            MAFFT alignment (reference first)
    MHCII_region_alignment.png/pdf alignment figure for the configured regions
    MHCII_identity_matrix.png/pdf  pairwise % sequence-identity heatmap
"""

import csv
import subprocess
import sys
from pathlib import Path

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import altair as alt
import seaborn as sns
import yaml
from Bio import SeqIO
from scipy import stats
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize

# Keep text as editable text (not outlined paths) in SVG output.
plt.rcParams["svg.fonttype"] = "none"

# Row order at the top of the figures; the first entry is also the residue
# numbering reference. Sequences not listed follow in input order. Names that
# are absent (e.g. because they were excluded) are simply ignored.
TOP_ORDER = ["tufted duck", "black swan", "bald eagle"]

# Clustal-style amino-acid coloring (by physicochemical class). Used as a
# fallback when no binding-escape data is supplied.
AA_COLORS = {
    "A": "#80b3e6", "I": "#80b3e6", "L": "#80b3e6", "M": "#80b3e6",
    "F": "#80b3e6", "W": "#80b3e6", "V": "#80b3e6", "C": "#80b3e6",  # hydrophobic
    "K": "#e6331a", "R": "#e6331a",                                   # positive
    "D": "#e600e6", "E": "#e600e6",                                   # negative
    "N": "#1acc1a", "Q": "#1acc1a", "S": "#1acc1a", "T": "#1acc1a",  # polar
    "G": "#e6994d",                                                   # glycine
    "P": "#cccc00",                                                   # proline
    "H": "#1ab3b3", "Y": "#1ab3b3",                                   # aromatic
}
GAP_COLOR = "#ffffff"      # gap (deletion relative to reference)
NODATA_COLOR = "#d9d9d9"   # residue/site not present in the binding data

# Diverging colormap for H5 HA binding escape: 0 (background) is white, positive
# escape is blue, negative (enhanced binding) is red.
ESCAPE_CMAP = plt.get_cmap("RdBu")


# --------------------------------------------------------------------------- #
# Small helpers
# --------------------------------------------------------------------------- #
def norm(name):
    """Normalize a sequence name for matching (lowercase, '_' -> space)."""
    return name.lower().replace("_", " ").strip()


def name_of(record):
    """Full sequence name (FASTA header text after '>', spaces preserved)."""
    return record.description


def parse_names(value):
    """Accept a string, comma-separated string, or list -> list of names."""
    if value is None:
        return []
    if isinstance(value, str):
        return [v.strip() for v in value.split(",") if v.strip()]
    return [str(v).strip() for v in value]


def parse_regions(spec):
    """Parse '31-39, 44-47' -> [(31, 39), (44, 47)] (1-based, inclusive)."""
    regions = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            start, end = chunk.split("-")
            regions.append((int(start), int(end)))
        else:
            regions.append((int(chunk), int(chunk)))
    return regions


def load_binding_escape(csv_path):
    """Read the DMS CSV -> (escape_lut, site_labels).

    `escape_lut` maps {(ORF_site, mutant): H5 HA binding escape}; ORF_site is the
    integer position in the reference numbering used for matching alignment
    columns. `site_labels` maps {ORF_site: site} where 'site' is the display
    label shown on the figure (e.g. "9(a)").

    Rows with a blank escape value are skipped. The reference (tufted duck)
    carries its wildtype residue at every site, for which the CSV stores an
    escape of 0, so reference squares come out at the colormap center.
    """
    lut = {}
    site_labels = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            orf_site = int(row["ORF_site"])
            site_labels[orf_site] = row["site"].strip()
            value = row["H5 HA binding escape"].strip()
            if value == "":
                continue
            lut[(orf_site, row["mutant"].strip())] = float(value)
    return lut, site_labels


def load_titers(csv_path, cells):
    """Read virus_titers.csv -> {virus (normalized): mean_RLUperuL} for `cells`.

    Only rows whose 'cells' column equals `cells` are kept.
    """
    titers = {}
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            if row["cells"].strip() != cells:
                continue
            value = row["mean_RLUperuL"].strip()
            if value == "":
                continue
            titers[norm(row["virus"])] = float(value)
    return titers


# --------------------------------------------------------------------------- #
# Sequence preparation and alignment
# --------------------------------------------------------------------------- #
def load_sequences(fasta, exclude_names, reference_name):
    """Read the FASTA, drop excluded sequences, and order rows for plotting.

    Returns the records ordered so TOP_ORDER names come first; the reference is
    forced to be the very first row. Raises if the reference was excluded or is
    missing.
    """
    records = list(SeqIO.parse(fasta, "fasta"))
    if not records:
        raise ValueError(f"No sequences found in {fasta}")

    excluded = {norm(n) for n in exclude_names}
    if norm(reference_name) in excluded:
        raise ValueError(f"Reference '{reference_name}' is in the exclude list.")

    kept = [r for r in records if norm(name_of(r)) not in excluded]
    dropped = [name_of(r) for r in records if norm(name_of(r)) in excluded]
    if dropped:
        print(f"Excluded: {', '.join(dropped)}")

    by_name = {norm(name_of(r)): r for r in kept}
    if norm(reference_name) not in by_name:
        raise ValueError(f"Reference '{reference_name}' not found in {fasta}")

    # Reference first, then the rest of TOP_ORDER, then everything else.
    order = [reference_name] + [n for n in TOP_ORDER if norm(n) != norm(reference_name)]
    front, seen = [], set()
    for n in order:
        if norm(n) in by_name:
            front.append(by_name[norm(n)])
            seen.add(norm(n))
    rest = [r for r in kept if norm(name_of(r)) not in seen]
    return front + rest


def mafft_align(records, out_fasta):
    """Align `records` with MAFFT, write to out_fasta, return aligned records.

    The aligned records are returned in the same row order as the input so the
    reference stays first.
    """
    out_path = Path(out_fasta)
    tmp_in = out_path.with_suffix(".input.fasta")
    SeqIO.write(records, tmp_in, "fasta")
    try:
        with open(out_path, "w") as out:
            subprocess.run(
                ["mafft", "--auto", "--anysymbol", str(tmp_in)],
                stdout=out,
                check=True,
            )
    finally:
        tmp_in.unlink(missing_ok=True)

    aligned = {norm(name_of(r)): r for r in SeqIO.parse(out_path, "fasta")}
    return [aligned[norm(name_of(r))] for r in records]


def ref_position_columns(ref_aligned_seq):
    """Map 1-based ungapped reference positions -> 0-based alignment columns."""
    mapping = {}
    res_num = 0
    for col, aa in enumerate(ref_aligned_seq):
        if aa != "-":
            res_num += 1
            mapping[res_num] = col
    return mapping


def summed_escape(aligned, regions, escape_lut):
    """Per-species summed H5 HA binding escape across the displayed region sites.

    Returns {sequence name: summed escape}, computed exactly like the per-row
    sums shown on the alignment figure. The reference (tufted duck) is wildtype
    at every site (escape 0), so its sum is 0; a residue with no matching escape
    entry contributes 0.
    """
    mapping = ref_position_columns(str(aligned[0].seq))
    col_to_pos = {col: pos for pos, col in mapping.items()}
    all_cols = []
    for start, end in regions:
        all_cols.extend(range(mapping[start], mapping[end] + 1))
    return {
        name_of(rec): sum(
            escape_lut.get((col_to_pos[col], str(rec.seq)[col]), 0.0)
            for col in all_cols if col in col_to_pos
        )
        for rec in aligned
    }


# --------------------------------------------------------------------------- #
# Figures
# --------------------------------------------------------------------------- #
def _text_color(rgba):
    """Black on light backgrounds, white on dark ones (luminance heuristic)."""
    lum = 0.299 * rgba[0] + 0.587 * rgba[1] + 0.114 * rgba[2]
    return "black" if lum > 0.55 else "white"


def plot_regions(aligned, regions, outfile_base, escape_lut=None, site_labels=None):
    """Plot the alignment restricted to the configured regions.

    Row 0 is the reference. Each region is a block of consecutive alignment
    columns spanning the reference residue range (so insertions in other
    sequences remain visible). Residues differing from the reference are
    outlined in black.

    If `escape_lut` ({(ORF_site, mutant): H5 HA binding escape}) is given, each
    square is colored by the escape value of the residue that sequence carries
    at that site (reference wildtype = 0, the colormap center). Otherwise
    squares use physicochemical amino-acid colors.

    Regions and matching use ORF_site (reference) numbering. If `site_labels`
    ({ORF_site: site}) is given, the per-column and region-header labels shown
    on the figure use the 'site' display label instead of the ORF_site number.
    """
    def site_label(orf_pos):
        """Display label for a reference (ORF) position."""
        if orf_pos is None:
            return None
        if site_labels and orf_pos in site_labels:
            return site_labels[orf_pos]
        return str(orf_pos)
    ref_aln = str(aligned[0].seq)
    mapping = ref_position_columns(ref_aln)
    labels = [name_of(r) for r in aligned]
    n_seqs = len(aligned)

    # Resolve each region to its alignment columns.
    blocks = []  # (start, end, [column indices])
    for start, end in regions:
        if start not in mapping or end not in mapping:
            raise ValueError(
                f"Region {start}-{end} exceeds reference length "
                f"({max(mapping)} residues)."
            )
        blocks.append((start, end, list(range(mapping[start], mapping[end] + 1))))

    col_to_pos = {col: pos for pos, col in mapping.items()}
    all_cols = [c for *_, cols in blocks for c in cols]

    # Symmetric color scale centered at 0, scaled to the escape values shown.
    escape_norm = None
    if escape_lut is not None:
        shown = [
            escape_lut[(col_to_pos[col], str(rec.seq)[col])]
            for col in all_cols if col in col_to_pos
            for rec in aligned
            if (col_to_pos[col], str(rec.seq)[col]) in escape_lut
        ]
        vmax = max((abs(v) for v in shown), default=1.0) or 1.0
        escape_norm = Normalize(vmin=-vmax, vmax=vmax)

    def cell_color(aa, ref_pos):
        """(facecolor, textcolor) for a residue `aa` at reference site `ref_pos`."""
        if aa == "-":
            return GAP_COLOR, "black"
        if escape_lut is not None:
            escape = escape_lut.get((ref_pos, aa)) if ref_pos is not None else None
            if escape is None:
                return NODATA_COLOR, "black"
            rgba = ESCAPE_CMAP(escape_norm(escape))
            return rgba, _text_color(rgba)
        return AA_COLORS.get(aa, "#dddddd"), "black"

    # Per-species sum of the escape values across the displayed sites (the
    # reference is wildtype everywhere, so its sum is 0).
    row_sums = None
    if escape_lut is not None:
        sums_by_name = summed_escape(aligned, regions, escape_lut)
        row_sums = [sums_by_name[name_of(rec)] for rec in aligned]

    gap = 1.0  # blank columns between region blocks
    sum_pad = 5.0 if row_sums is not None else 1.0  # right-hand room for the sums
    total_cols = sum(len(cols) for *_, cols in blocks)
    fig_w = (max(len(l) for l in labels) * 0.11
             + (total_cols + gap * len(blocks)) * 0.32 + sum_pad * 0.32 + 0.5)
    fig_h = (n_seqs + 2.0) * 0.32 + 0.6

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.set_xlim(0, total_cols + gap * len(blocks) + sum_pad)
    ax.set_ylim(0, n_seqs + 2)
    ax.invert_yaxis()
    ax.axis("off")

    x = 0
    for start, end, cols in blocks:
        ax.text(x + len(cols) / 2.0, -0.4,
                f"{site_label(start)}–{site_label(end)}",
                ha="center", va="bottom", fontsize=9, fontweight="bold")

        for ci, col in enumerate(cols):
            ref_aa = ref_aln[col]
            ref_pos = col_to_pos.get(col)
            if ref_pos is not None:
                ax.text(x + ci + 0.5, 0.55, site_label(ref_pos),
                        ha="center", va="bottom", fontsize=6, rotation=90,
                        color="#444444")

            for row, rec in enumerate(aligned):
                aa = str(rec.seq)[col]
                differs = row != 0 and aa != ref_aa
                facecolor, textcolor = cell_color(aa, ref_pos)
                ax.add_patch(mpatches.Rectangle(
                    (x + ci, row + 1), 1, 1,
                    facecolor=facecolor,
                    edgecolor="black" if differs else "#cccccc",
                    linewidth=1.4 if differs else 0.3,
                ))
                if aa != "-":
                    ax.text(x + ci + 0.5, row + 1.5, aa, ha="center", va="center",
                            fontsize=7, color=textcolor,
                            fontweight="bold" if differs else "normal")

        x += len(cols) + gap

    # Per-species summed escape, in a column to the right of the last block.
    if row_sums is not None:
        sx = x + 0.2
        ax.text(sx, -0.4, "Σ escape\n(shown sites)", ha="left", va="bottom",
                fontsize=8, fontweight="bold")
        for row, total in enumerate(row_sums):
            ax.text(sx, row + 1.5, f"{total:+.2f}", ha="left", va="center",
                    fontsize=8, fontweight="bold" if row == 0 else "normal")

    for row, label in enumerate(labels):
        ax.text(-0.4, row + 1.5, label, ha="right", va="center", fontsize=8,
                fontweight="bold" if row == 0 else "normal")

    ax.set_title(
        f"MHCII alpha chain alignment — binding-site regions "
        f"(numbering: {labels[0]})",
        fontsize=11, pad=20,
    )

    if escape_norm is not None:
        sm = ScalarMappable(norm=escape_norm, cmap=ESCAPE_CMAP)
        cbar = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
        cbar.set_label("H5 HA binding escape", fontsize=9)

    _save(fig, outfile_base)


def plot_titer_correlation(aligned, regions, escape_lut, titers, outfile_base,
                           cells_label):
    """Altair scatter of per-species summed DMS escape vs measured titer.

    x = summed H5 HA binding escape across the displayed region sites
    y = mean_RLUperuL in `cells_label` (plotted on a log scale, since the titers
        span several orders of magnitude)

    Only species present in both the alignment and the titer table are shown.
    Pearson r is computed on (escape, log10 titer) and shown in the subtitle.
    """
    sums_by_name = summed_escape(aligned, regions, escape_lut)
    points = [
        (name, escape, titers[norm(name)])
        for name, escape in sums_by_name.items()
        if norm(name) in titers
    ]
    missing = [name for name in sums_by_name if norm(name) not in titers]
    if missing:
        print(f"No {cells_label} titer for: {', '.join(missing)} (omitted)")
    if len(points) < 2:
        print("Fewer than 2 species with both escape and titer; "
              "skipping correlation plot.")
        return

    df = pd.DataFrame(points, columns=["species", "escape", "titer"])
    escape = df["escape"].to_numpy()
    log_titer = np.log10(df["titer"].to_numpy())
    pearson_r, pearson_p = stats.pearsonr(escape, log_titer)

    # Least-squares fit on the plotted (escape, log10 titer) scale.
    slope, intercept = np.polyfit(escape, log_titer, 1)

    # Pad the domains so no point sits on an axis (y padded in log space).
    x_pad = 0.12 * (np.ptp(escape) or 1.0)
    x_domain = [escape.min() - x_pad, escape.max() + x_pad]
    y_pad = 0.12 * (np.ptp(log_titer) or 1.0)
    y_domain = [10 ** (log_titer.min() - y_pad), 10 ** (log_titer.max() + y_pad)]

    xs = np.linspace(*x_domain, 100)
    fit_df = pd.DataFrame({"escape": xs, "titer": 10 ** (slope * xs + intercept)})

    x_enc = alt.X(
        "escape:Q",
        title="Σ H5 HA binding escape (shown sites)",
        scale=alt.Scale(zero=False, nice=False, domain=x_domain),
    )
    y_enc = alt.Y(
        "titer:Q",
        title=f"mean RLU/µL ({cells_label})",
        scale=alt.Scale(type="log", nice=False, domain=y_domain, clamp=True),
    )

    fit = alt.Chart(fit_df).mark_line(
        color="#9e9e9e", strokeDash=[6, 4], strokeWidth=1.5,
    ).encode(x=x_enc, y=alt.Y("titer:Q"))

    pts = alt.Chart(df).mark_point(
        size=140, filled=True, opacity=0.9, color="#0081A7",
        stroke="black", strokeWidth=0.7,
    ).encode(
        x=x_enc, y=y_enc,
        tooltip=["species:N", alt.Tooltip("escape:Q", format=".3f"),
                 alt.Tooltip("titer:Q", format=",.0f")],
    )

    labels = alt.Chart(df).mark_text(
        align="left", dx=8, dy=-6, fontSize=12,
    ).encode(x=x_enc, y=y_enc, text="species:N")

    chart = (fit + pts + labels).properties(
        width=340, height=260,
        title=alt.TitleParams(
            "Summed DMS escape vs measured titer",
            subtitle=f"Pearson r = {pearson_r:.2f}  (p = {pearson_p:.2g}, "
                     f"n = {len(df)})",
            fontSize=14, subtitleFontSize=11, anchor="start",
        ),
    ).configure_view(stroke=None).configure_axis(
        labelFontSize=10, titleFontSize=11, grid=False,
    )

    exts = ("html", "png", "svg")
    for ext in exts:
        kw = {"ppi": 200} if ext == "png" else {}
        chart.save(f"{outfile_base}.{ext}", **kw)
    print(f"Wrote {', '.join(f'{outfile_base}.{e}' for e in exts)}")


def plot_identity_matrix(aligned, outfile_base):
    """Lower-triangle heatmap of pairwise % identity.

    Styled to match mhcii_sequence_identity_heatmap.py: a white -> #0081A7
    colormap, only the lower triangle shown, white bold annotations, a 50-100
    color range and thin cell borders.
    """
    labels = [name_of(r) for r in aligned]
    seqs = [str(r.seq) for r in aligned]
    n = len(aligned)

    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            mat[i, j] = _percent_identity(seqs[i], seqs[j])
    identity_df = pd.DataFrame(mat, index=labels, columns=labels)

    fig, ax = plt.subplots(figsize=(5, 4.5))

    # White -> #0081A7 colormap, matching the reference figure.
    custom_cmap = LinearSegmentedColormap.from_list(
        "custom_blue", ["#FFFFFF", "#0081A7"], N=256
    )

    # Show only the lower triangle (mask the strict upper triangle).
    mask = np.triu(np.ones_like(identity_df.values, dtype=bool), k=1)

    sns.heatmap(
        identity_df,
        mask=mask,
        annot=True,
        fmt=".1f",
        cmap=custom_cmap,
        vmin=30,
        vmax=100,
        square=True,
        linewidths=0.2,
        cbar_kws={"label": "Sequence Identity (%)", "shrink": 0.7},
        annot_kws={"size": 10, "color": "white", "weight": "bold"},
        ax=ax,
    )

    ax.set_title("MHCII Alpha Chain Sequence Identity",
                 fontsize=12, fontweight="bold", pad=10)
    ax.set_xlabel("MHCII Sequences", fontsize=8, fontweight="bold")
    ax.set_ylabel("MHCII Sequences", fontsize=8, fontweight="bold")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)
    plt.setp(ax.get_yticklabels(), rotation=0, fontsize=7)

    _save(fig, outfile_base)


def _percent_identity(seq_a, seq_b):
    """% identity over alignment columns where neither sequence is a gap."""
    matches = compared = 0
    for x, y in zip(seq_a, seq_b):
        if x == "-" or y == "-":
            continue
        compared += 1
        matches += x == y
    return 100.0 * matches / compared if compared else 0.0


def _save(fig, outfile_base):
    fig.tight_layout()
    exts = ("png", "pdf", "svg")
    for ext in exts:
        fig.savefig(f"{outfile_base}.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {', '.join(f'{outfile_base}.{e}' for e in exts)}")


# --------------------------------------------------------------------------- #
def main():
    config_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("config.yaml")
    config = yaml.safe_load(config_path.read_text())

    fasta = Path(config["MHCII_sequences"])
    reference_name = config["compared_to"]
    regions = parse_regions(config["Tufted_duck_MHCII_binding_sites"])
    exclude_names = parse_names(config.get("exclude"))

    escape_lut = None
    site_labels = None
    binding_csv = config.get("Tufted_duck_MHCII_binding_data")
    if binding_csv:
        escape_lut, site_labels = load_binding_escape(Path(binding_csv))
        print(f"Loaded H5 HA binding escape for {len(escape_lut)} (ORF_site, mutant) pairs")

    print(f"Reference (top row, numbering): {reference_name}")
    print(f"Regions (reference numbering):  {regions}")

    records = load_sequences(fasta, exclude_names, reference_name)
    aligned = mafft_align(records, "MHCII_aligned.fasta")
    print(f"Aligned {len(aligned)} sequences; "
          f"reference length {max(ref_position_columns(str(aligned[0].seq)))} residues")

    plot_regions(aligned, regions, "MHCII_region_alignment",
                 escape_lut=escape_lut, site_labels=site_labels)
    plot_identity_matrix(aligned, "MHCII_identity_matrix")

    titer_csv = config.get("virus_titer_data")
    if titer_csv and escape_lut is not None:
        cells = config.get("titer_cells", "M3-Wigeon-H5")
        titers = load_titers(Path(titer_csv), cells)
        print(f"Loaded {cells} titers for {len(titers)} viruses")
        plot_titer_correlation(aligned, regions, escape_lut, titers,
                               "MHCII_escape_vs_titer", cells)


if __name__ == "__main__":
    main()
