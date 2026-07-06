---
aside: false
---

# Effects of mutations in tufted duck MHC-II on binding to H5 HA

This page shows how mutations to tufted duck MHC-II affect binding to H5 HA, as assessed by the ability of soluble HA trimer to neutralize MHC-II pseudovirus infectivity.The H5 HA used in these experiments comes from A/American Wigeon/South Carolina/USDA-000345-001/2021 (H5N1) strain.

[[toc]]

## Heatmap of how mutations affect H5 HA binding

In the interactive chart below, the line plot shows the total effects of mutations at each site in tufted duck MHC-II on H5 HA binding, and the heatmap shows the effects of individual mutations.

The heatmap shows the effects of individual mutations: positive values (blue) indicate mutations increase binding (as assessed by increased HA neutralization of pseudovirus) and negative values (orange) indicate mutations decrease binding.
Light gray in the heatmap means a mutation was not measured, and dark gray means the mutation was so deleterious for MHC-II-mediated cell entry that it was not possible to reliably estimate its effect on HA binding.
The `x` for each site indicates the wildtype amino acid in the tufted duck MHC-II. Note that the numbering in these plots refers to sites in tufted duck MHC-II ectodomain (with `(a)` indicating alpha chain and `(b)` indicating beta chain ectodomain numbering. Note that these site numbers are not equivalent to homologous site numbers in human MHC-II)
Mouse over mutations on the heatmap for details.

Below the chart are interactive options to adjust parameters on the chart.
Click the box in the upper right to expand the chart to full page.

<Figure caption="Effects of mutations in tufted duck MHC-II on H5 HA binding">
    <Altair :showShadow="true" :spec-url="'htmls/H5_HA_binding_mut_effect.html'"></Altair>
</Figure>

Additional relevant links:

 - [standalone link to the chart shown above](htmls/H5_HA_binding_mut_effect.html){target="_self"}
 - [chart showing effects of mutations in tufted duck MHC-II on entry in H5 HA expressing cells](htmls/entry_all_cells_overlaid.html){target="_self"}

## MHC-II structure colored by mutation effects on HA binding

Below is a structure of tufted duck MHC-II colored by the total effect of all mutations at each site on H5 HA binding. A peptide has been modeled in from a different structure (`PDB:8JRJ`) for visual orientation purposes is shown in green.
Red indicates sites where most mutations decrease HA binding and blue indicates sites where most mutations increase HA binding

<video autoplay muted loop playsinline width="100%" style="display: block;">
  <source src="./MHCII.mp4" type="video/mp4">
</video>

## Numerical values of mutation effects on MHC-II binding
For numerical data, see the following CSVs:
  - [Effects of mutations on H5 HA binding after recommended QC](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/summaries/HA_binding.csv) 
  - [Measured effects on MHC-II binding with full QC details but QC not pre-applied to numerical values](https://github.com/dms-vep/Tufted-duck-MHCII-DMS/blob/master/results/HA_binding/averages/HA_binding_mut_effect.csv) (only use this file if you understand QC filters)
