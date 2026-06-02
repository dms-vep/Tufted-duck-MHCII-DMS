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