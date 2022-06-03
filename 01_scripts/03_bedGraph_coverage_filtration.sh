#!/bin/bash
# bedGraph_coverage_filtration.sh
#srun -c 1 --mem 10G -p medium --time 7-00:00:00 -J 03_cov_filt -o 03_cov_filt_%j.log ./01_scripts/03_bedGraph_coverage_filtration.sh &

INPUT="04_masked_bedGraphs"
OUTPUT="05_quick_filtered_bedGraphs"

for i in $(ls -1 "$INPUT"/*dedup_CpG_merged.bedGraph.gz | cut -d "." -f1) 
do
    echo $(basename $i)
    gunzip -c "$INPUT"/"$(basename $i)".dedup_CpG_merged.bedGraph.gz | 
		awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $3, $4, $5, $6, ($5+$6)}' |
		awk -v FS="\t" -v OFS="\t" '{ if (4<$7 && 101>$7) {print} } ' |
    gzip -c - > "$OUTPUT"/"$(basename $i)""_quickfiltered.bedGraph.gz"
done