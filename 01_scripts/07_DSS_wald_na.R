#!/usr/bin/env Rscript
#conda activate R411 (newer version of DSS)
#srun -c 5 --mem 150G -p small --time 1-00:00:00 -J 07_DSS -o 07_DSS_wald_%j.log Rscript ./01_scripts/07_DSS_wald_na.R &
# script for two-group comparison with Wald test
# install_github("haowulab/DSS")
# conda install -c conda-forge r-base=4.1.0

### WALD MODEL SPECIFICATIONS
#metadata="99_sample_info/hugo_capelin_sample_metadata_europe.txt"
#output_prefix="DSS_results_europe"

#metadata="99_sample_info/hugo_capelin_sample_metadata_northamerica.txt"
#output_prefix="DSS_results_northamerica"

metadata="99_sample_info/hugo_capelin_sample_metadata_northamerica_noBB65.txt"
output_prefix="DSS_results_northamerica_noBB65"

#metadata="99_sample_info/hugo_capelin_sample_metadata_northamerica_noBB65_noRIG1_noMAK2_DSEA_15.txt"
#output_prefix="DSS_results_northamerica_noBB65_noRIG1_noMAK2_DSEA_15"

groups=3	# number of column of metadata file with test factor
p_val=0.01
delta=0.1	# minimum percent difference in methylation between groups
min_len=200
min_CGs=10
merge_distance=50
pct_sig_CGs=0.5
smooth=TRUE	# TRUE or FALSE - true is recommended. FALSE for debugging (faster)
nCPUs=5

library(DSS)
library(data.table)
library(dplyr)

# load files
path = "07_DSS_input/"
files <- read.table(metadata,header=T)[,1]
make_df = lapply(paste0(path, files, ".txt.gz"), function (x) fread(x, header=T))

# make BSseq object
BSobj = makeBSseqData(make_df, basename(files))
BSobj

print(paste("Saving .RData file."))
save.image(paste0("08_DSS_results/", output_prefix, ".RData"))
#load(paste0("08_DSS_results/", output_prefix, ".RData"))


# experimental design - create lists of sample names for groups to plug into Wald test
sample_info <- read.table(metadata, header=T)
unique(sample_info[groups])[[1]]
group1 <- sample_info[1][sample_info[groups]==unique(sample_info[groups])[[1]][1]]
group2 <- sample_info[1][sample_info[groups]==unique(sample_info[groups])[[1]][2]]


# DML test - parallelized across nCPUs
print(paste("DML testing. Comparing", unique(sample_info[groups])[[1]][1], "and", unique(sample_info[groups])[[1]][2]))
system.time(DMLfit <- DMLtest(BSobj, group1=group1, group2=group2, smoothing=smooth, nCPUs))
#fwrite(DMLfit, file="DMLfit.txt", sep="\t", quote=FALSE)
#head(DMLfit)
DML_test = callDML(DMLfit, delta=delta, p.threshold=p_val)


# filter DMLs by p-value and write tables
print(paste("Filtering DMLs."))
sig_DMLs <- DML_test[which(DML_test$fdr < p_val),]
fwrite(sig_DMLs, file=paste0("08_DSS_results/", output_prefix, "_", colnames(sample_info[groups]), "_DMLs_FDR", p_val, "_delta", delta, ".txt"), sep="\t", quote=FALSE, row.names=FALSE)


print(paste("Calling DMRs."))
DMR_test = callDMR(DMLfit, delta=delta, p.threshold=p_val, minlen=min_len, minCG=min_CGs, dis.merge=merge_distance, pct.sig=pct_sig_CGs)
write.table(DMR_test, file=paste0("08_DSS_results/", output_prefix, "_", colnames(sample_info[groups]), "_DMRs_p", p_val, "_minlen", min_len, "_", min_CGs, "minCGs_", merge_distance, "bpmerge_", pct_sig_CGs, "pctsigCGs.txt"), sep="\t", quote=FALSE, row.names=FALSE)


print(paste("Getting heatmap beta values."))
DMRegions <- getMeth(BSobj, DMR_test, type="raw", what="perRegion")
betas <- cbind(DMR_test[,1:3], DMRegions)
rownames(betas)<-NULL
betas <- betas %>% mutate_all(~ifelse(is.nan(.), NA, .))
write.table(betas, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[2], "_DMRs_p", p_val, "_minlen", min_len, "_", min_CGs, "minCGs_", merge_distance, "bpmerge_", pct_sig_CGs, "pctsigCGs_heatmapbetas.txt"), sep="\t", quote=FALSE, row.names=FALSE)
