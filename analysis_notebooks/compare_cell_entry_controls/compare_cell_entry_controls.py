"""Correlate functional-effect columns across the four cell-entry-control spreadsheets.

Reads the four CSV paths from ``config.yaml``, merges them on a concatenated
``site`` + ``wildtype`` + ``mutant`` key, and builds an interactive Altair figure
of pairwise scatterplots between the ``effect`` columns. Each scatterplot shows a
Pearson correlation coefficient that updates with the sliders, and the sliders
filter on the ``times_seen`` / ``n_selections`` columns of the
``VSVG_infection_multi_mut`` spreadsheet.

Run with the ``dms-vep-pipeline-3`` conda environment, e.g.::

    conda run -n dms-vep-pipeline-3 python compare_cell_entry_controls.py
"""

from pathlib import Path

import altair as alt
import pandas as pd
import yaml

# the merge key column derived from the four spreadsheets
KEY = "site_wildtype_mutant"

# spreadsheet config keys -> short column suffix used in the merged dataframe
SHEETS = {
    "VSVG_infection_single_mut": "inf_single",
    "VSVG_RT_single_mut": "RT_single",
    "VSVG_infection_multi_mut": "inf_multi",
    "VSVG_RT_multi_mut": "RT_multi",
}

# human-readable axis titles per spreadsheet
TITLES = {
    "inf_single": "infection single-mut effect",
    "RT_single": "RT single-mut effect",
    "inf_multi": "infection multi-mut effect",
    "RT_multi": "RT multi-mut effect",
}

# the spreadsheet whose times_seen / n_selections drive the sliders
SLIDER_SHEET = "inf_multi"

# pairwise comparisons between the four effect columns
PAIRS = [
    ("inf_single", "RT_single"),
    ("inf_single", "inf_multi"),
    ("inf_single", "RT_multi"),
    ("RT_single", "inf_multi"),
    ("RT_single", "RT_multi"),
    ("inf_multi", "RT_multi"),
]


def resolve_path(relpath, config_dir):
    """Resolve a config-relative path, trying a few sensible base directories."""
    for base in (config_dir, config_dir.parent, config_dir.parent.parent):
        candidate = (base / relpath).resolve()
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"could not locate {relpath!r} relative to {config_dir}")


def load_and_merge(config_path):
    """Load the four spreadsheets and outer-merge them on the concatenated key."""
    config_dir = config_path.parent
    with open(config_path) as f:
        config = yaml.safe_load(f)

    merged = None
    for sheet_key, suffix in SHEETS.items():
        csv_path = resolve_path(config[sheet_key], config_dir)
        df = pd.read_csv(csv_path)
        df[KEY] = df["wildtype"].astype(str) + df["site"].astype(str) + df["mutant"].astype(str)

        cols = {KEY: KEY, "effect": f"effect_{suffix}"}
        # carry times_seen / n_selections only from the slider spreadsheet
        if suffix == SLIDER_SHEET:
            cols["times_seen"] = "times_seen_slider"
            cols["n_selections"] = "n_selections_slider"

        df = df[list(cols)].rename(columns=cols)
        merged = df if merged is None else merged.merge(df, on=KEY, how="outer")

    return merged


def make_scatter(merged, xsuffix, ysuffix, filt):
    """One scatterplot of two effect columns with a reactive correlation label."""
    xcol, ycol = f"effect_{xsuffix}", f"effect_{ysuffix}"
    xtitle, ytitle = TITLES[xsuffix], TITLES[ysuffix]

    base = alt.Chart(merged).transform_filter(filt)

    points = base.mark_circle(size=18, opacity=0.4, color="#1f77b4").encode(
        x=alt.X(f"{xcol}:Q", title=xtitle, scale=alt.Scale(padding=15)),
        y=alt.Y(f"{ycol}:Q", title=ytitle, scale=alt.Scale(padding=15)),
        tooltip=[
            alt.Tooltip(f"{KEY}:N", title="wildtype/site/mutant"),
            alt.Tooltip(f"{xcol}:Q", title=xtitle, format=".3f"),
            alt.Tooltip(f"{ycol}:Q", title=ytitle, format=".3f"),
            alt.Tooltip("times_seen_slider:Q", title="times_seen (inf multi)"),
            alt.Tooltip("n_selections_slider:Q", title="n_selections (inf multi)"),
        ],
    )

    # Pearson r computed on the filtered rows so it tracks the sliders
    corr = (
        base.transform_joinaggregate(mx=f"mean({xcol})", my=f"mean({ycol})")
        .transform_calculate(
            dxy=f"(datum['{xcol}'] - datum.mx) * (datum['{ycol}'] - datum.my)",
            dx2=f"(datum['{xcol}'] - datum.mx) * (datum['{xcol}'] - datum.mx)",
            dy2=f"(datum['{ycol}'] - datum.my) * (datum['{ycol}'] - datum.my)",
        )
        .transform_aggregate(sxy="sum(dxy)", sx2="sum(dx2)", sy2="sum(dy2)", n="count()")
        .transform_calculate(r="datum.sxy / sqrt(datum.sx2 * datum.sy2)")
        .transform_calculate(
            label="'r = ' + format(datum.r, '.2f') + '  (n = ' + datum.n + ')'"
        )
        .mark_text(align="left", baseline="top", fontSize=12, fontWeight="bold", color="black")
        .encode(x=alt.value(6), y=alt.value(6), text=alt.Text("label:N"))
    )

    return (points + corr).properties(width=240, height=240)


def main():
    config_path = Path(__file__).resolve().parent / "config.yaml"
    merged = load_and_merge(config_path)

    # data is larger than Altair's default 5000-row cap
    alt.data_transformers.disable_max_rows()

    # sliders that filter on the inf-multi spreadsheet's times_seen / n_selections
    times_seen_param = alt.param(
        name="times_seen_min",
        value=0,
        bind=alt.binding_range(min=0, max=50, step=0.5, name="min times_seen (inf multi): "),
    )
    n_selections_param = alt.param(
        name="n_selections_min",
        value=1,
        bind=alt.binding_range(min=1, max=2, step=1, name="min n_selections (inf multi): "),
    )

    def filt_for(xsuffix, ysuffix):
        xcol, ycol = f"effect_{xsuffix}", f"effect_{ysuffix}"
        return (
            f"isValid(datum['{xcol}']) && isValid(datum['{ycol}'])"
            " && datum.times_seen_slider >= times_seen_min"
            " && datum.n_selections_slider >= n_selections_min"
        )

    charts = [make_scatter(merged, x, y, filt_for(x, y)) for x, y in PAIRS]

    figure = (
        alt.concat(*charts, columns=2)
        .add_params(times_seen_param, n_selections_param)
        .properties(
            title="Correlation of functional effects across cell-entry controls"
        )
        .configure_title(anchor="start")
    )

    out_path = config_path.parent / "compare_cell_entry_controls.html"
    figure.save(out_path)
    print(f"merged {len(merged)} mutations; wrote {out_path}")


if __name__ == "__main__":
    main()
