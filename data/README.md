# Input data
This subdirectory contains input data used by the pipeline.

## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps several different numbering shcemes for tufted duck MHC-II. Columns in the spreadsheet include: *sequential_site* (sequential numbering of MHC-II 1,2,3...), *reference_site* (alpha (a), beta (b) or T2A linker (T2A) domains numbered individually. Negative values indicate signal peptide sequence and positive values are mature protein domains numbered sequentially), *structure_site* (site numbers matching numbering in tufted duck alpha fold MHC-II structure, here undercores a or b mean site is in alpha or beta chain, respectively and for each chain sites are numbered starting first position after the signal peptide, sites not part of the extodomain are labeled as NA). *ORF_site* (numbered starting position 22, which is the first amino acid in the ectodomain in alpha chain. This numbering is useful when converting reference_site to position in the open reading frame of a plasmid construct).

## Mutation-type classification
[data/mutation_design_classification.csv](data/mutation_design_classification.csv) classifies mutations into the different categories of designed mutations.
Has columns *sequential_site*, *amino_acid*, and *mutation_type*.

## Neutralization standard barcodes
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) barcodes for the neutralization standards.
Has columns *barcode* and *name*, giving the barcode and name of this neutralization standard set.

## Barcode runs
[barcode_runs.csv](barcode_runs.csv) contains all samples and paths to sequencing files. It has the following format:

 - `sample`: sample name
 - `library`: name of library
 - `date`: date of sequencing
 - `fastq_R1`: path to one more FASTQ R1 sequencing files, multiple files should be semicolon-delimited

## Configuration for analyzing functional effects of mutations
[func_effects_config.yml](func_effects_config.yml) has the configuration for analyzing functional effects of mutations.
The format is explained within the file.

