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
#		awk -v FS="\t" -v OFS="\t" '{ if (4<$7) {print} } ' |
#    gzip -c - > "$OUTPUT"/"$(basename $i)""_SNPsexcluded2_noscaffolds_coveragefiltered.bedGraph.gz"
    gzip -c - > "$OUTPUT"/"$(basename $i)""_quickfiltered.bedGraph.gz"
done




#awk -v FS="\t" -v OFS="\t" 'NR {print $1, $2, $3, $4, $5, $6, ($5+$6)}' CD17_SNPsexcluded2_noscaffolds_whitelisted.bedGraph > CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum.bedGraph
#awk -v FS="\t" -v OFS="\t" '{ if (4<$7 && $7<101) {print} } '  CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum.bedGraph > CD17_SNPsexcluded2_noscaffolds_whitelisted_testsum_covfilter.bedGraph
