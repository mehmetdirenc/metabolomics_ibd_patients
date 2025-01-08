library(MetaboAnalystR)
# library(MSnbase)
# library(xcms)
# library(dplyr)
# library(tibble)
# library(pheatmap)
# library(DESeq2)
# library(limma)

data <- read.csv("/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", sep = ",", header = TRUE, row.names = 1)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, "/home/direnc/inputs/mahana/ms_ms/metabolomics_statistics/patient_abx_muc2_pre_post_serum.csv", "colu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-SanityCheckData(mSet)
mSet<-FilterVariable(mSet, "F", 20, "iqr", 5, "mean", 0)
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "SumNorm", "LogNorm", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-FC.Anal(mSet, 2.0, 0, FALSE)
mSet<-PlotFC(mSet, "fc_0_", "png", 72, width=NA)
# mSet$analSet$fc$fc.log
# view(mSet$analSet$fc$fc.log)
mSet<-Ttests.Anal(mSet, F, 0.05, FALSE, TRUE, "none", FALSE)
mSet<-PlotTT(mSet, "tt_0_", "png", 72, width=NA)
mSet<-PlotTT(mSet, imgName = "tt_0_", format = "png", dpi = 72, width=NA)
ttest_results <- mSet$analSet$tt$sig.mat
view(ttest_results)


# x <- rownames(data)
# group <- as.factor(as.character(data["Label", ]))
# data <- data[-1, ]  # Remove the Label row for processing
# x <- rownames(data)
# data <- as.data.frame(lapply(data, as.numeric))  # Ensure numeric values
# rownames(data) <- x  # Assign metabolite names
#
#
# mSet <- InitDataObjects("spec", "stat", FALSE)
# mSet <- Read.TextData(mSet, data, grp = group, format = "rowu")
#
#
# mSet <- Normalization(mSet, "NULL", "LogNorm", "NULL", ratio = FALSE)
#
# mSet <- PerformFoldChange(mSet, 1.5, 0)
#
# mSet <- Ttests.Anal(mSet)
#
# mSet <- AdjustPvals(mSet, "fdr")
#
# fc_table <- mSet$analSet$fc$resTable
# ttest_table <- mSet$analSet$tt$resTable

#
# mSet <- InitDataObjects("stat", "pktable", FALSE)
#
# # Load the data
# mSet <- Read.TextData(mSet, "/home/direnc/inputs/mahana/ms_ms/test.csv", "rowu", "disc")
#
# # Check your data
# # View(mSet$dataSet$orig)
#
# mSet <- FilterVariable(mSet, "iqr")
#
# # Normalize data (e.g., log transformation, Pareto scaling)
# mSet <- PreparePrenormData(mSet)
# mSet <- Normalization(mSet, "NULL", "LogNorm", "ParetoNorm", ratio = FALSE, ref = NULL)



# print("asd")
# mzml_files <- list.files(path = "/home/direnc/inputs/mahana/ms_ms/metaboanalystr_test/", pattern = "*.mzML", full.names = TRUE)
# raw_data <- readMSData(files = mzml_files, mode = "onDisk")
# sample_groups <- c("blank", "blank", "cal6pre", "cal6pre", "cal6pre")
# group_param <- PeakDensityParam(sampleGroups = sample_groups, bw = 30, minFraction = 0.5, binSize = 0.005)
#
# peaks <- findChromPeaks(raw_data, param = CentWaveParam())
#
# peaks <- groupChromPeaks(peaks, param = group_param, msLevel = 1L)
# grouped_data <- chromPeaks(peaks)
# head(grouped_data)  # Preview grouped peaks
#
# feature_table <- featureValues(peaks, value = "into")  # Use 'into' for integrated intensity
# feature_table <- as.data.frame(feature_table)
#
# annotation_file_path <- "/home/direnc/inputs/mahana/ms_ms/metaboanalystr_test/compounds.csv"
# metabolite_annotations <- read.table(annotation_file_path, header = FALSE, sep = ",")
# colnames(metabolite_annotations) <- c("mz", "compound_name")
#
# # Step 6: Prepare feature_table for joining by adding row names as feature_id
# feature_table <- feature_table %>%
#     rownames_to_column(var = "feature_id") %>%
#     left_join(metabolite_annotations, by = "feature_id")
#
#
# # Step 8: View the updated feature_table
# head(feature_table)
