#!/usr/bin/env Rscript
#~/syidan/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/qc/Human_ThyroidEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR_hichip_qc.tsv
#~/syidan/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/qc/Human_Th17_unknown_WT_unknown_HiCHIP_standard_maps_specific_qc.tsv
#/mnt/citadel2/research/syidan/Data/Processed/HiCHIP_GSE/ATAC/merged.snakePipes.out/deepTools_qc/estimateReadFiltering/Human_Treg_unknown_WT_unknown_ATAC_standard_filtering_estimation.txt

# Load required libraries
library(dplyr)
library(ggplot2)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Please provide exactly 5 arguments: input_dir, filename_pattern, metric_name, output_dir, output_name_prefix\nUsage: Rscript plot_column_barplot.R /path/to/input/dir *pattern metric_name /path/to/output/dir prefix")
}

# Define arguments
input_dir <- args[1]
file_pattern <- args[2]
metric_name <- args[3]
output_dir <- args[4]
output_name_prefix <- args[5]

# Validate input and output directories
if (!dir.exists(input_dir)) {
  stop(paste("Input directory does not exist:", input_dir))
}
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created output directory:", output_dir, "\n")
}

# Define output files
output_txt <- file.path(output_dir, paste0("summary_", output_name_prefix, ".txt"))
output_plot <- file.path(output_dir, paste0("summary_", output_name_prefix, "_barplot.pdf"))

# Find all files matching the pattern
tsv_files <- list.files(path = input_dir, pattern = file_pattern, full.names = TRUE)
if (length(tsv_files) == 0) {
  stop(paste("No files matching", file_pattern, "found in", input_dir))
}

# Read and extract specified metric and sample name from each file
if (file_pattern == "flagstat.txt"){
  data_list <- lapply(tsv_files, function(file) {
    print(file)
    df <- first_line <- readLines(file, n = 1)
    #print(df)
    sample_name <- basename(file)
    metric_value <- as.numeric(gsub(" .*$", "", df))
    # Split sample name by underscore and take 2nd and 5th parts
    sample_parts <- unlist(strsplit(as.character(sample_name), "_"))
    if (length(sample_parts) < 4) {
      warning(paste("Sample name in", file, "has fewer than 5 parts:", sample_name))
      short_name <- sample_name
    } else {
      short_name <- gsub("_unknown","",paste(sample_parts[2], sample_parts[4], sample_parts[5], sep = "_"))
    }
    return(data.frame(Sample = short_name, Metric_Value = as.numeric(gsub("[\\(\\)%]", "", metric_value)), Original_Sample = sample_name))
  })
} else {
  data_list <- lapply(tsv_files, function(file) {
    print(file)
    df <- read.table(file, header = FALSE, sep = "\t", fill = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    print(df)
    #if (length(grep(metric_name, df[[1]])) != 0) {
    #  stop(paste(metric_name, "metric not found in", file))
    #}
    sample_name <- basename(file)
    metric_value <- df[grep(metric_name, df[[1]]), 2]
    # Split sample name by underscore and take 2nd and 5th parts
    sample_parts <- unlist(strsplit(as.character(sample_name), "_"))
    if (length(sample_parts) < 4) {
      warning(paste("Sample name in", file, "has fewer than 5 parts:", sample_name))
      short_name <- sample_name
    } else {
      short_name <- gsub("_unknown","",paste(sample_parts[2], sample_parts[4], sample_parts[5], sep = "_"))
    }
    return(data.frame(Sample = short_name, Metric_Value = as.numeric(gsub("[\\(\\)%]", "", metric_value)), Original_Sample = sample_name))
  })
}
# Combine into a single data frame
data <- bind_rows(data_list) %>% rename(!!metric_name := Metric_Value)

# Write summary to text file
write.table(data, output_txt, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote summary to", output_txt, "\n")

# Create barplot
p <- ggplot(data, aes(x = reorder(Sample, .data[[metric_name]]), y = .data[[metric_name]])) +
  geom_col(fill = "#FFA4A4") +
  labs(title = "",
       x = "Sample",
       y = metric_name) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Add horizontal line at 20M for Valid_Pairs
if (metric_name == "Valid_Pairs") {
  p <- p + geom_hline(yintercept = 20000000, linetype = "dashed", color = "black")
} else if (metric_name == "FRiP") {
  p <- p + geom_hline(yintercept = 20, linetype = "dashed", color = "black")
}

# Save plot
ggsave(output_plot, plot = p, width = 15, height = 8, dpi = 300)
cat("Wrote barplot to", output_plot, "\n")

cat("Processing completed successfully\n")