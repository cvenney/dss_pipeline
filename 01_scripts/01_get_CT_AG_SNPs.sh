#!/bin/bash
# srun -c 1 --mem 1G -p small --time 1-00:00:00 -J 01_CTpoly -o 01_CTpoly_%j.log ./01_scripts/01_get_CT_AG_SNPs.sh &

INPUT="02_SNPs/sites_all_maf0.01_pctind0.2"

awk '($3 == "C" && $4 == "T") || ($3 == "T" && $4 == "C") || ($3 == "G" && $4 == "A") || ($3 == "A" && $4 == "G"){
    	print $1 "\t" $2 - 1 "\t" $2 "\t" $1 "_" $2 "\t" $3 "\t" $4
    }' "$INPUT" > "$INPUT"_CT_AG_snps.bed