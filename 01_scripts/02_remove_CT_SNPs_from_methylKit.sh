#!/bin/bash
# srun -c 4 --mem 10G -p small --time 1-00:00:00 -J 02_kitSNP -o 02_kitSNP_%j.log ./01_scripts/02_remove_CT_SNPs_from_methylKit.sh &

INPUT="/project/lbernatchez/users/clven7/hugo_capelin/bwa-meth_pipeline/07_methyl_dackel/"
OUTPUT="04_masked_bedGraphs"
NCPUS=4

module load bedtools

# Negative intersect with bedtools
for file in $(ls "$INPUT"/*.methylKit.gz | perl -pe 's/\.methylKit\.gz//g')
do
    name=$(basename $file)
    
    echo "Filtering sample: $file"    
    
    gunzip -c "$INPUT"/"$name".methylKit.gz | 
    awk 'BEGIN{OFS="\t"}(NR!=1){print $2, $3 -1, $3, $1, $4, $5, $6, $7}' |
    bedtools intersect -a stdin -b "/project/lbernatchez/users/clven7/hugo_capelin/dss_pipeline/02_SNPs/sites_all_maf0.01_pctind0.2_CT_AG_snps.bed" -wa -v |
    awk 'BEGIN{OFS="\t"; print "chrBase", "chr", "base", "strand", "coverage", "freqC", "freqT"} {print $4, $1, $3, $5, $6, $7, $8}' |
    gzip > "$OUTPUT"/"$name".methylKit.gz

done
