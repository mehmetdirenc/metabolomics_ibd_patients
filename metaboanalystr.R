library(MetaboAnalystR)
# library(MSnbase)
# library(xcms)
# library(dplyr)
# library(tibble)
# library(pheatmap)
# library(DESeq2)
# library(limma)

# data <- read.csv("/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", sep = ",", header = TRUE, row.names = 1)
input_dir <- "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_ours"
output_dir <- "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics_ours/"

heatmap_data <- data.frame()

file_list <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
ttest_results_list <- list()
for (file_path in file_list) {
  # Extract file name without extension for output naming
  file_name <- tools::file_path_sans_ext(basename(file_path))

  # Initialize MetaboAnalystR data objects
  mSet <- InitDataObjects("pktable", "stat", FALSE)
  mSet <- Read.TextData(mSet, file_path, "colu", "disc")
  mSet <- SanityCheckData(mSet)
  mSet <- ReplaceMin(mSet)
  mSet <- SanityCheckData(mSet)
  mSet <- FilterVariable(mSet, "F", 20, "iqr")
  mSet <- PreparePrenormData(mSet)
  mSet <- Normalization(mSet, "SumNorm", "LogNorm", "AutoNorm", ratio = FALSE, ratioNum = 20)
  # mSet <- PlotNormSummary(mSet, paste0(output_dir, "norm_", file_name, "_"), "png", 72)
  # mSet <- PlotSampleNormSummary(mSet, paste0(output_dir, "snorm_", file_name, "_"), "png", 72)
  #### FOLD CHANGE DIRECTION IS A/B E.G. CTR_IBD3_CAECUM CITRULLINE = 0.03 A HAS LOWER LEVELS THAN B ####
  #### FOLD CHANGE DIRECTION IS A/B E.G. CTR_IBD3_CAECUM D-GLUCOSAMINE 6-PHOSPHATE = 21.2 A HAS HIGHER LEVELS THAN B ####
  mSet <- FC.Anal(mSet, 2, 0, FALSE)
  # mSet <- PlotFC(mSet, paste0(output_dir, "fc_", file_name, "_"), "png", 72)

  # Perform T-tests
  mSet <- Ttests.Anal(mSet, nonpar = F, threshp = 0.05, paired = FALSE, equal.var = TRUE, "fdr", FALSE)
  # mSet <- PlotTT(mSet, imgName = paste0(output_dir, "tt_", file_name, "_"), format = "png", dpi = 72)
  ttest_results <- mSet$analSet$tt$sig.mat
  compound_names <- rownames(ttest_results)
  ttest_results <- cbind(Compounds = compound_names, ttest_results)
  rownames(ttest_results) <- NULL
  # colnames(ttest_results) <- c("Compounds", "t.stat", "p.value", "-log10(p)", "FDR")
  all_fold_changes <- mSet$analSet$fc$fc.all
  fold_changes <- all_fold_changes
  ttest_results <- as.data.frame(ttest_results)
  ttest_results$FoldChange <- fold_changes[match(ttest_results$Compounds, names(fold_changes))]
  # compound_names_fc <- rownames(fold_changes)
  # fold_changes <- cbind(Compounds = compound_names_fc, fold_changes)
  # ttest_results$Fold_Change <- fold_changes[rownames(ttest_results)]

  # Save individual T-test results with fold changes
  # write.csv(ttest_results, file = paste0(output_dir, "ttest_results_", file_name, ".csv"))

  # Combine results for the heatmap
  ttest_results_list[[file_name]] <- ttest_results
  output_file <- paste0(output_dir, file_name, "_ttest_results.tsv")
  write.table(ttest_results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # write.csv(heatmap_data, file = paste0(output_dir, "combined_ttest_results.csv"))
  # heatmap_data <- rbind(heatmap_data, ttest_results)
}

# Save the combined results
# write.csv(heatmap_data, file = paste0(output_dir, "combined_ttest_results.csv"))

# mSet<-InitDataObjects("pktable", "stat", FALSE)
# mSet<-Read.TextData(mSet, "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", "colu", "disc")
# mSet<-SanityCheckData(mSet)
# mSet<-ReplaceMin(mSet)
# mSet<-SanityCheckData(mSet)
# mSet<-FilterVariable(mSet, "F", 20, "iqr")
# mSet<-PreparePrenormData(mSet)
# mSet<-Normalization(mSet, "SumNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
# mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
# mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
# mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
# mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
# mSet$analSet$fc$fc.log
# # view(mSet$analSet$fc$fc.log)
# # mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "none", FALSE)
# mSet<-Ttests.Anal(mSet, nonpar=T, threshp=0.05, paired=FALSE, equal.var=FALSE, "none", FALSE)
# mSet<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
# mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)
# ttest_results <- mSet$analSet$tt$sig.mat
# View(ttest_results)




# data <- read.csv("/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", sep = ",", header = TRUE, row.names = 1)
#
# mSet<-InitDataObjects("pktable", "stat", FALSE)
# mSet<-Read.TextData(mSet, "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", "colu", "disc")
# mSet<-SanityCheckData(mSet)
# mSet<-ReplaceMin(mSet)
# mSet<-SanityCheckData(mSet)
# mSet<-FilterVariable(mSet, "F", 20, "iqr")
# mSet<-PreparePrenormData(mSet)
# mSet<-Normalization(mSet, "SumNorm", "LogNorm", "AutoNorm", ratio=FALSE, ratioNum=20)
# mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
# mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
# mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
# mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
# mSet$analSet$fc$fc.log
# # view(mSet$analSet$fc$fc.log)
# # mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "none", FALSE)
# mSet<-Ttests.Anal(mSet, nonpar=T, threshp=0.05, paired=FALSE, equal.var=FALSE, "none", FALSE)
# mSet<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
# mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)
# ttest_results <- mSet$analSet$tt$sig.mat
# View(ttest_results)
