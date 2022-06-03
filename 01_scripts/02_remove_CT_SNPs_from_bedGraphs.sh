#!/bin/bash
# srun -c 1 --mem 2G -p small --time 1-00:00:00 -J 02_bedSNP -o 02_bedSNP_%j.log ./01_scripts/02_remove_CT_SNPs_from_bedGraphs.sh &

INPUT="/project/lbernatchez/users/clven7/hugo_capelin/bwa-meth_pipeline/07_methyl_dackel/"
OUTPUT="04_masked_bedGraphs"
NCPUS=4

module load bedtools

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.bedGraph.gz | perl -pe 's/\.bedGraph\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".bedGraph.gz | 
    bedtools intersect -a stdin -b "/project/lbernatchez/users/clven7/hugo_capelin/dss_pipeline/02_SNPs/sites_all_maf0.01_pctind0.2_CT_AG_snps.bed" -wa -v | 
    gzip > "$OUTPUT"/"$name".bedGraph.gz

done
