# filter_DMRs_by_areastat.R

nDMRs=2000
dmrfile="08_DSS_results/DSS_results_europe_habitat_DMRs_p0.01_minlen200_10minCGs_50bpmerge_0.5pctsigCGs.txt"

library(dplyr)
library(tools)

dmr <- read.table(dmrfile, header=T)

dmrslice <- dmr %>% slice_max(abs(areaStat), n = nDMRs)

write.table(dmrslice, paste0("08_DSS_results/DMR_subsets/", file_path_sans_ext(basename(dmrfile)), "_top", nDMRs, "DMRs.txt"), quote=F)
