#!/bin/bash
#bedgraph_coverage_info.sh
# get mean, median, and SD coverage for each sample
# srun -o  bedgraph_coverage_info_%j.log ./01_scripts/bedgraph_coverage_info.sh &

INPUT=06_fully_filtered_bedGraphs

for i in $(ls -1 "$INPUT"/*_fully_filtered.bedGraph.gz | cut -d "_" -f4) 
do
	echo $(basename $i)
	gunzip -c "$INPUT"/"$(basename $i)"_*.gz |
	awk '{sum+=$7} END { print sum/NR }'
	
done