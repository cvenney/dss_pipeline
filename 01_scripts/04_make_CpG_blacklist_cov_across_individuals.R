#!/usr/bin/env Rscript
#conda activate R403b
#usage: Rscript ./01_scripts/make_CpG_blacklist_cov_across_individuals.R 4 sample_info/cliff_dwarf.txt sample_info/cliff_normal.txt 05_coverage_filtered_bedGraphs/03_make_CpG_blacklist_cov_across_individuals.R
#srun -c 1 --mem 40G -p medium --time 7-00:00:00 -J 04_make_blacklist -o 04_make_blacklist_%j.log Rscript ./01_scripts/04_make_CpG_blacklist_cov_across_individuals.R 11 99_sample_info/hugo_capelin_sample_IDs.txt 05_quick_filtered_bedGraphs/blacklist_max_11_missing.bed &

library(data.table)
library(dplyr)

argv <- commandArgs(T)
max_missing_ind <- as.numeric(argv[1])
sample_file <- argv[2]
output <- argv[3]

max_missing_ind <- 11
sample_file <- "99_sample_info/hugo_capelin_sample_IDs.txt"
output <- "blacklist_max_11_missing.txt"


vector_ind<-read.table(sample_file,header=F)[,1]
#read ind 1
bedgraph<-fread(paste0("05_quick_filtered_bedGraphs/", vector_ind[1],"_quickfiltered.bedGraph.gz"))[,1:4]
colnames(bedgraph)<-c("chr","start","stop",vector_ind[1])

#loop from ind 2 to ind n to get full matrix
for (i in 2 : length(vector_ind))
{
print(vector_ind[i])
#bedgraph_i<-fread(paste0(vector_ind[i],"_SNPsexcluded2_noscaffolds_coveragefiltered.bedGraph.gz")
bedgraph_i<-fread(paste0("05_quick_filtered_bedGraphs/", vector_ind[i],"_quickfiltered.bedGraph.gz"))
bedgraph_4<-bedgraph_i[,1:4]
colnames(bedgraph_4)<-c("chr","start","stop",vector_ind[i])
bedgraph<-full_join(bedgraph,bedgraph_4)
}

count_na <- function(x) sum(is.na(x))

bedgraph_filtered <- bedgraph %>%
	mutate(count_na = apply(., 1, count_na)) %>%
	filter(count_na > max_missing_ind)
dim(bedgraph)
dim(bedgraph_filtered)

print(paste("unfiltered CpGs:", nrow(bedgraph)-1))
print(paste("CpGs in blacklist (insufficient coverage across individuals):", nrow(bedgraph_filtered)-1))
print(paste("total CpGs after filtering", ((nrow(bedgraph)-1) - (nrow(bedgraph_filtered)-1))))

write.table(bedgraph_filtered[,1:3], output, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
