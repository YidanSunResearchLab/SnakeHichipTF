#!/usr/bin/env Rscript

library(data.table)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript hichip_qc.R <sample> <output_file> <mapping_stats> <interaction_stats> <peaks> [optional_inputs...]")
}

# Fixed arguments (first 5)
sample_name <- args[1]
output_file <- args[2]
mapping_stats <- args[3]  # hicpro/.../sample_hg38.bwt2pairs.pairstat
interaction_stats <- args[4]  # hicpro/.../sample_allValidPairs.mergestat
peaks_file <- args[5]

# Optional arguments (remaining args, 0 to 4 files)
optional_inputs <- args[6:length(args)]  # Everything after the 5th argument
hicdcplus_loops <- NA
hichipper_loops <- NA
fithichip_loops <- NA
maps_loops <- NA

# Assign optional inputs based on file patterns
if (length(optional_inputs) > 0) {
  for (i in seq_along(optional_inputs)) {
    if (grepl("hicdcplus", optional_inputs[i])) hicdcplus_loops <- optional_inputs[i]
    else if (grepl("hichipper", optional_inputs[i])) hichipper_loops <- optional_inputs[i]
    else if (grepl("fithichip", optional_inputs[i])) fithichip_loops <- optional_inputs[i]
    else if (grepl("maps", optional_inputs[i])) maps_loops <- optional_inputs[i]
  }
}

# 1. Sequencing and Mapping Stats (from .pairstat)
mapping <- fread(mapping_stats, sep = "\t", header = FALSE)
total_reads <- as.integer(mapping[V1 == "Total_pairs_processed", V2])
mapped_reads <- as.integer(mapping[V1 == "Unique_paired_alignments", V2]) + 
                as.integer(mapping[V1 == "Multiple_pairs_alignments", V2])
valid_pairs_count <- as.integer(mapping[V1 == "Reported_pairs", V2])

# 2. Interaction Stats (from .mergestat)
interactions <- fread(interaction_stats, sep = "\t", header = FALSE)
valid_interaction_rmdup <- as.integer(interactions[V1 == "valid_interaction_rmdup", V2])
cis_interactions <- as.integer(interactions[V1 == "cis_interaction", V2])
trans_interactions <- as.integer(interactions[V1 == "trans_interaction", V2])
cis_short_range <- as.integer(interactions[V1 == "cis_shortRange", V2])
cis_long_range <- as.integer(interactions[V1 == "cis_longRange", V2])
cis_pct <- (cis_interactions / valid_interaction_rmdup) * 100
trans_pct <- (trans_interactions / valid_interaction_rmdup) * 100

# 3. Peak Stats
peaks <- fread(peaks_file)
num_peaks <- nrow(peaks)
peak_widths <- peaks$V3 - peaks$V2
median_width <- median(peak_widths)
#mean_signal <- mean(-log10(peaks$V7))  # Assuming column 7 is p-value

# 4. Interaction Metrics (handle missing files with NA)
# HICDCPLUS
hicdcplus_loops_count <- NA
hicdcplus_median_dist <- NA
if (!is.na(hicdcplus_loops)) {
  hicdcplus <- fread(hicdcplus_loops, header = TRUE)
  hicdcplus_loops_count <- nrow(hicdcplus)
  hicdcplus_median_dist <- median(abs(as.numeric(hicdcplus$startI) - as.numeric(hicdcplus$startJ)))
}

# HICHIPPER
hichipper_loops_count <- NA
hichipper_median_dist <- NA
if (!is.na(hichipper_loops)) {
  hichipper <- fread(hichipper_loops, header = FALSE)
  hichipper_loops_count <- nrow(hichipper)
  hichipper_median_dist <- median(abs(as.numeric(hichipper$V2) - as.numeric(hichipper$V5)))
}

# FITHICHIP
fithichip_loops_count <- NA
fithichip_median_dist <- NA
if (!is.na(fithichip_loops)) {
  fithichip <- fread(fithichip_loops, header = TRUE)
  #print(head(fithichip))
  fithichip_loops_count <- nrow(fithichip)
  fithichip_median_dist <- median(abs(fithichip$s1 - fithichip$s2))
}

# MAPS
maps_loops_count <- NA
maps_median_dist <- NA
if (!is.na(maps_loops)) {
  maps <- fread(maps_loops, header = TRUE)
  #print(head(maps))
  maps_loops_count <- nrow(maps)
  maps_median_dist <- median(abs(maps$start1 - maps$start2))  # Adjust if bin1_mid/bin2_mid
}

# 5. Overlap (simplified example, only for available tools)
common_loops <- NA
peak_overlap_pct <- NA
if (sum(!is.na(c(hicdcplus_loops, hichipper_loops, fithichip_loops, maps_loops))) > 1) {
  loop_sets <- list()
  if (!is.na(hicdcplus_loops)) loop_sets[["hicdcplus"]] <- fread(hicdcplus_loops)$V1
  if (!is.na(hichipper_loops)) loop_sets[["hichipper"]] <- fread(hichipper_loops)$V1
  if (!is.na(fithichip_loops)) loop_sets[["fithichip"]] <- fread(fithichip_loops)$V1
  if (!is.na(maps_loops)) loop_sets[["maps"]] <- fread(maps_loops)$V1
  common_loops <- length(Reduce(intersect, loop_sets))

#  overlap_checks <- c()
#  if (!is.na(maps_loops)) overlap_checks <- c(overlap_checks, any(fread(maps_loops)$V2 %in% peaks$V2 | fread(maps_loops)[,5] %in% peaks$V2))
#  if (!is.na(hichipper_loops)) overlap_checks <- c(overlap_checks, any(fread(hichipper_loops)$V1 %in% peaks$V2 | fread(hichipper_loops)$V4 %in% peaks$V2))
# if (!is.na(fithichip_loops)) overlap_checks <- c(overlap_checks, any(fread(fithichip_loops)$V2 %in% peaks$V2 | fread(fithichip_loops)[,5] %in% peaks$V2))
#  if (!is.na(hicdcplus_loops)) overlap_checks <- c(overlap_checks, any(fread(hicdcplus_loops)$V2 %in% peaks$V2 | fread(hicdcplus_loops)[,5] %in% peaks$V2))
#  peak_overlap_pct <- mean(overlap_checks) * 100
}

print(paste("hicdcplus_loops:", hicdcplus_loops))
print(paste("hichipper_loops:", hichipper_loops))
print(paste("fithichip_loops:", fithichip_loops))
print(paste("maps_loops:", maps_loops))

# 6. Write QC File
qc <- data.frame(
  Sample = sample_name,
  Total_Reads = total_reads,
  Mapped_Reads = mapped_reads,
  Valid_Pairs = valid_pairs_count,
  Valid_Interactions_Rmdup = valid_interaction_rmdup,
  Cis_Pct = cis_pct,
  Trans_Pct = trans_pct,
  Cis_Short_Range = cis_short_range,
  Cis_Long_Range = cis_long_range,
  Num_Peaks = num_peaks,
  Peak_Median_Width = median_width,
  #Peak_Mean_Signal = mean_signal,
  MAPS_Loops = maps_loops_count,
  MAPS_Median_Distance = maps_median_dist,
  HICHIPPER_Loops = hichipper_loops_count,
  HICHIPPER_Median_Distance = hichipper_median_dist,
  FITHICHIP_Loops = fithichip_loops_count,
  FITHICHIP_Median_Distance = fithichip_median_dist,
  HICDCPLUS_Loops = hicdcplus_loops_count,
  HICDCPLUS_Median_Distance = hicdcplus_median_dist
  #Common_Loops = common_loops,
  #Peak_Overlap_Pct = peak_overlap_pct
)
fwrite(t(qc), output_file, sep = "\t", quote = FALSE, col.names = FALSE, row.names = TRUE)