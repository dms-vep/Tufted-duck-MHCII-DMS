# Deep mutational scanning of tufted duck MHC-II for mutation effects on H5 HA binding

See [Dadonaite et al (2026)](https://doi.org/10.64898/2026.07.17.738765) for details on this study.

This repository has the code and data for analysis of tufted duck MHC-II deep mutational scanning, including how mutations affect binding to a H5 HA. The tufted duck MHC-II sequence used comed from Genscript accessions XP_032061117.1 and XP_032061025.1 for alpha and beta chain, respectively.

For rendering of key results and an easy-to-interpret summary, see the documentation of the analysis at [https://dms-vep.org/Tufted-duck-MHCII-DMS/](https://dms-vep.org/Tufted-duck-MHCII-DMS/).

The key results file with the processed measurements is [results/summaries/HA_binding.csv](results/summaries/HA_binding.csv).

## Organization of this repo

### `dms-vep-pipeline-3` submodule

Most of the analysis is done by the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3), which was added as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) to this pipeline via:

    git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This added the file [.gitmodules](.gitmodules) and the submodule [dms-vep-pipeline-3](dms-vep-pipeline-3), which was then committed to the repo.
Note that if you want a specific commit or tag of [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) or to update to a new commit, follow the [steps here](https://stackoverflow.com/a/10916398), basically:

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory, and add and commit the updated `dms-vep-pipeline-3` submodule.
You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) that you commit back to that repo.

### Additional custom rules
Additional custom rules outside the standard pipeline are defined in [custom_rules.smk](custom_rules.smk). Notebooks associated with these custom rules are in [./pipeline_notebooks/][pipeline_notebooks] subforlder. 

### Configuration and running the pipeline
The configuration for the pipeline is in [config.yaml](config.yaml) and the files in [./data/](data) referenced therein.
To run the pipeline, do:

    snakemake -j 8 --software-deployment-method conda -s dms-vep-pipeline-3/Snakefile

To run on the Hutch cluster via [slurm](https://slurm.schedmd.com/), you can run the file [run_Hutch_cluster.bash](run_Hutch_cluster.bash):

    sbatch -c 8 run_Hutch_cluster.bash

### Input data
Input data for the pipeline are in [./data/](data). 

### Results and documentation
The results of running the pipeline are placed in [./results/](results).
Only some of these results are tracked to save space (see [.gitignore](.gitignore)).

The pipeline builds HTML documentation for the pipeline in `./results/docs` and `./results/publish_docs`.
To visualize these docs via GitHub Pages, run:

    dms-vep-pipeline-3/publish_docs_gh-pages.sh

This pushes the docs to the *gh-pages* branch, we can be viewed on GitHub Pages at [https://dms-vep.org/Tufted-duck-MHCII-DMS/](https://dms-vep.org/Tufted-duck-MHCII-DMS/).

### Non-pipeline analyses
Analysis notebooks and scripts that are not part of the main dms-vep pipeline are in [./analysis_notebooks/](analysis_notebooks) subfolder.

