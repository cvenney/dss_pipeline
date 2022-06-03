#!/usr/bin/env Rscript
#conda activate R403b
#srun -c 1 --mem 100G -p medium --time 7-00:00:00 -J 07_DSS -o 07_DSS_%j.log Rscript ./01_scripts/07_DSS_glm.R &
# currenty written for 2 main terms and one interaction term!

### GLM MODEL SPECIFICATIONS
metadata="99_sample_info/hugo_capelin_sample_metadata.txt"
model_formula=as.formula("~habitat + continent + habitat:continent")
output_prefix="DSS_results"
p_val=0.001
min_len=200
min_CGs=10
merge_distance=50
pct_sig_CGs=0.5


library(DSS)
library(data.table)

# load files
path = "07_DSS_input/"
files = list.files(path=path, pattern=paste0("*.txt.gz"))
make_df = lapply(paste0(path, files), function (x) fread(x, header=T))

# make BSseq object
#BSobj=makeBSseqData(dat=make_df, sampleNames=files)
BSobj = makeBSseqData(make_df, basename(files))
BSobj

print(paste("Saving .RData file."))
save.image(paste0("08_DSS_results/", output_prefix, ".RData"))
#load(paste0("08_DSS_results/", output_prefix, ".RData"))


# experimental design
sample_info <- read.table(metadata, header=T)
vars <- sample_info[,1:3]


# DML test
# colnames(DMLfit$X)[1] is the intercept. Add terms as needed ([5] and beyond)
print(paste("DML testing."))
DMLfit = DMLfit.multiFactor(BSobj, design=vars, formula=model_formula)
colnames(DMLfit$X)
one <- colnames(DMLfit$X)[2]
two <- colnames(DMLfit$X)[3]
three <- colnames(DMLfit$X)[4]


DMLtest.t1 = DMLtest.multiFactor(DMLfit, coef=one)
DMLtest.t2 = DMLtest.multiFactor(DMLfit, coef=two)
DMLtest.t3 = DMLtest.multiFactor(DMLfit, coef=three)

# if you want to sort DMLs
#ix=sort(DMLtest.ecotype[,"pvals"], index.return=TRUE)$ix
#head(DMLtest.ecotype[ix,])

# if you want a list of all CpG sites analyzed and data on them
#write.table(DMLtest.ecotype, file="08_DSS_results/cliff_DMLs_all_CpGs.txt", sep="\t", quote=FALSE)

# filter DMLs by p-value and write tables
print(paste("Filtering DMLs."))
sig_DMLs_1 <- DMLtest.t1[which(DMLtest.t1$fdrs < p_val),]
fwrite(sig_DMLs_1, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[2], "_DMLs_FDR", p_val, ".txt"), sep="\t", quote=FALSE)

sig_DMLs_2 <- DMLtest.t2[which(DMLtest.t2$fdrs < p_val),]
fwrite(sig_DMLs_2, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[3], "_DMLs_FDR", p_val, ".txt"), sep="\t", quote=FALSE)

sig_DMLs_3 <- DMLtest.t3[which(DMLtest.t3$fdrs < p_val),]
fwrite(sig_DMLs_3, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[4], "_DMLs_FDR", p_val, ".txt"), sep="\t", quote=FALSE)


print(paste("Calling DMRs."))
DMRtest_1 = callDMR(DMLtest.t1, p.threshold=p_val, minlen=min_len, minCG=min_CGs, dis.merge=merge_distance, pct.sig=pct_sig_CGs)
write.table(DMRtest_1, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[2], "_DMRs_p", p_val, "_minlen", min_len, "_", min_CGs, "minCGs_", merge_distance, "bpmerge_", pct_sig_CGs, "pctsigCGs.txt"), sep="\t", quote=FALSE)

DMRtest_2 = callDMR(DMLtest.t2, p.threshold=p_val, minlen=min_len, minCG=min_CGs, dis.merge=merge_distance, pct.sig=pct_sig_CGs)
write.table(DMRtest_2, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[3], "_DMRs_p", p_val, "_", min_CGs, "minCGs_", merge_distance, "bpmerge_", pct_sig_CGs, "pctsigCGs.txt"), sep="\t", quote=FALSE)

DMRtest_3 = callDMR(DMLtest.t3, p.threshold=p_val, minlen=min_len, minCG=min_CGs, dis.merge=merge_distance, pct.sig=pct_sig_CGs)
write.table(DMRtest_3, file=paste0("08_DSS_results/", output_prefix, "_", colnames(DMLfit$X)[4], "_DMRs_p", p_val, "_", min_CGs, "minCGs_", merge_distance, "bpmerge_", pct_sig_CGs, "pctsigCGs.txt"), sep="\t", quote=FALSE)
