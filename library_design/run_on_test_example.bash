#!/bin/bash

python gga_codon_muts_oligo_design.py \
    --tiles_csv tufted_duck_MHCII_assembly_fragments.csv \
    --mutations_to_make_csv mutations_to_make.csv \
    --output_oligos_fasta GGA_oPool.fa \
    --avoid_motifs 'CGTCTC'
