"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""

rule entry_binding_correlations:
    """Compare entry and binding."""
    input:
        nb="pipeline_notebooks/compare_entry_and_binding.ipynb",
        summary_df="results/summaries/HA_binding.csv",
    output:
        nb="results/pipeline_notebooks/compare_entry_and_binding.ipynb",
        corr_chart="results/compare_cell_entry/corr_chart.html",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    log:
        "results/logs/entry_binding_correlations.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p summary_df {input.summary_df} \
            -p corr_chart {output.corr_chart} \
            &> {log}
        """

docs["Additional plots"] = {
    "Correlation between binding and entry":
        rules.entry_binding_correlations.output.corr_chart   
}

# Make row-wrapped heatmaps -------------------------------------------------------------

# read configuration for wrapped heatmaps
with open("data/wrapped_heatmap_config.yaml") as f:
    wrapped_heatmap_config = yaml.YAML(typ="safe", pure=True).load(f)


rule wrapped_heatmap:
    """Make row-wrapped heatmaps."""
    input:
        data_csv=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]["data_csv"],
    output:
        chart_html="results/wrapped_heatmaps/{wrapped_hm}_wrapped_heatmap.html",
    params:
        params_dict=lambda wc: wrapped_heatmap_config[wc.wrapped_hm]
    log:
        notebook="results/pipeline_notebooks/wrapped_heatmap_{wrapped_hm}.ipynb",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml"),
    notebook:
        "pipeline_notebooks/wrapped_heatmap.py.ipynb"

docs["Row-wrapped heatmaps"] = {
    "Heatmap HTMLs" : {
        wrapped_hm: rules.wrapped_heatmap.output.chart_html.format(wrapped_hm=wrapped_hm)
        for wrapped_hm in wrapped_heatmap_config
    }
}