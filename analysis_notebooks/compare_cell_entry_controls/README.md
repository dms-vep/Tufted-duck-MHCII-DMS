# Compare cell-entry controls

This code plots the correlation between cell-entry effects measured using
different VSVG controls.

The two controls differ in how the VSVG sample was prepared for Illumina library
prep:

- **RT** control: a VSVG sample made by extracting RNA from a rescued VSVG stock
  and reverse transcribing it before doing Illumina library prep.
- **infection** control: made by infecting cells with rescued VSVG virus and
  extracting non-integrated viral genomes using a miniprep kit before performing
  Illumina library prep.

## Usage

Paths to the four functional-effect spreadsheets (single- and multi-mutant
effects for each of the two controls) are set in `config.yaml`. Run with the
project conda environment:

```bash
conda run -n dms-vep-pipeline-3 python compare_cell_entry_controls.py
```

This writes `compare_cell_entry_controls.html`, an interactive Altair figure of
pairwise scatterplots between the `effect` columns of the four spreadsheets. The
spreadsheets are merged on a concatenated `wildtype` + `site` + `mutant` key,
which is shown on hover. Each scatterplot displays a Pearson correlation
coefficient, and two sliders filter the points on the `times_seen` and
`n_selections` columns of the `VSVG_infection_multi_mut` spreadsheet.
