---
layout: home

hero:
  name: "Effects of mutations to a tufted duck MHC-II on its interaction with H5 HA"
  tagline: "Pseudovirus deep mutational scanning of how mutations to tufted duck MHC-II affects its interaction with H5 HA"
  image: logo_MHCII.png
features:
  - title: Binding to H5 HA
    details: Effects of tufted duck MHC-II mutations on binding to H5 HA
    link: /binding
  - title: Effects of H5 HA mutations on MHC-II interaction
    details: Deep mutational scanning of how mutations to H5 HA affect its interaction with tufted duck MHC-II
    link: https://dms-vep.org/Flu-H5N1-American-Wigeon-2021-HA-tufted-duck-MHCII-DMS/
---

## Overview

This website provides interactive visualizations and links to numerical data from pseudovirus deep mutational scanning measuring how mutations to a tufted duck MHC-II affect its interaction with H5 HA. The tufted duck MHC-II used in these experiments has  Genbank accession numbers `XP_032061117.1` and `XP_032061025.1` for alpha and beta chains, respectively. The H5 HA used in these experiments comes from A/American Wigeon/South Carolina/USDA-000345-001/2021 (H5N1) strain.

For details about the study, see [Dadonaite et al. (2026)]() [**ADD CITATION**].

Visualizations and data can be accessed by clicking the gray boxes above for each measurement, namely:
 - [Effects of tufted duck mutations on binding to H5 HA](binding)
 - [Effects of H5 HA mutations on MHC-II interaction](https://dms-vep.org/Flu-H5N1-American-Wigeon-2021-HA-tufted-duck-MHCII-DMS/)

For numerical values of the effects of tufted duck MHC-II mutations on interaction with HA, see [this CSV](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/summaries/HA_binding.csv).

For the full computer code and all numerical data files, see [the GitHub repository](https://github.com/dms-vep/Tufted-duck-MHCII-DMS).
For full documentation of the computational pipeline, see the [Appendix](appendix.html){target="_self"}.

## Inverted pseudotyping
Note that these measurements are made using inverted pseudotyping of lentiviral particles, as described in [Dadonaite et al. (2026)]() [**ADD CITATION**].

The [numerical results]((https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/summaries/HA_binding.csv) therefore include both the effects of MHC-II mutations on H5 HA binding (what you are likely interested in) as well as on inverted pseudovirus entry (this is a proxy for MHC-II folding, with mutations have negative cell entry effects if they greatly decrease cell entry).

## MHC-II sequence numbering
The tufted duck MHC-II protein sequence used in these experiments (Genbank accession numbers `XP_032061117.1` and `XP_032061025.1` for alpha and beta chains, respectively) is [here](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/gene_sequence/protein.fasta), and involves the alpha and beta chains connected by a 2A linker.
The numbering used here is 1, 2, ... numbering of the **ectodomains** (see [here](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/data/site_numbering_map.csv) for how that numbering scheme relates to sequential numbering of the [2A-linked protein](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/gene_sequence/protein.fasta); note that this is **not** the same as the standard numbering scheme used for human MHC-II proteins.

## Biosafety
These experiments use lentiviral particles that inverse pseudotyped with MHC-II.
These lentiviral particles encode no viral proteins, and therefore are not pathogens capable of causing disease.

See the biosafety statement in [Dadonaite et al. (2026)]() [**ADD CITATION**] for more details.
