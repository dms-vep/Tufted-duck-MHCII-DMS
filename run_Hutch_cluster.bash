#!/bin/bash
#
#SBATCH -c 4

snakemake -j 4 --software-deployment-method conda -s dms-vep-pipeline-3/Snakefile
