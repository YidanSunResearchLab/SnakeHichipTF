#!/usr/bin/env Rscript

# Load required libraries
library(diffloop)
library(ggplot2)

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check for required arguments (expecting 7 arguments)
if (length(args) != 7) {
  stop("Usage: Rscript hichipper_diffloop.R <control_name> <test_name> <control_files> <test_files> <input_dir> <output_dir> <fdr>")
}

# Assign arguments to variables
control_name <- args[1]
test_name <- args[2]
control_files <- strsplit(args[3], " ")[[1]]
test_files <- strsplit(args[4], " ")[[1]]
input_dir <- args[5]
output_dir <- args[6]
fdr <- as.numeric(args[7])

# Print parsed options for debugging
print(paste("Control group:", control_name))
print(paste("Test group:", test_name))

# Create output directory if it doesn’t exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Step 1: Load all BEDPE files into a single loops object
sample_names <- c(control_files, test_files)
class(sample_names)
print(control_files)
print(test_files)
print(sample_names)

# Use loopsMake to load loops (assuming BEDPE format from HiChIP data)
message("Executing: all_loops <- loopsMake(beddir = ", input_dir, ", samples = c(", 
        paste(shQuote(sample_names, type = "cmd"), collapse = ", "), "), mergegap = 0)")
all_loops <- loopsMake(beddir = input_dir, samples = paste0(sample_names, ".filt.intra"))


# Step 2: Assign group labels to distinguish control vs. test conditions
group_labels <- c(rep("control", length(control_files)), rep("test", length(test_files)))
all_loops <- updateLDGroups(all_loops, group_labels)
tmp=summary(all_loops)
head(tail(tmp))
dim(tmp)

# Step 3: Quality Control (QC) Steps from the Vignette
# 3.1 Remove self-ligated loops (loopWidth < 5000 bp)
all_loops <- subsetLoops(all_loops, all_loops@rowData$loopWidth >= 5000)

# 3.2 Apply Mango correction for Hi-C bias (FDR as specified)
all_loops_corrected <- mangoCorrection(all_loops, FDR = fdr)

# Step 4: Generate PCA plot
pca_plot <- pcaPlot(all_loops_corrected) + 
  ggtitle("PCA of HiChIP Samples") + 
  theme_minimal()

# Save PCA plot to the output directory
pca_file <- paste0(output_dir, "/diffloop_pca_plot.tiff")
ggsave(pca_file, plot = pca_plot, width = 4, height = 3, dpi = 300)
message("PCA plot saved to ", pca_file)

# Step 5: Perform differential analysis with quickAssoc
message("Running differential analysis between ", control_name, " and ", test_name)
diff_result <- quickAssoc(all_loops_corrected)

# Step 6: Extract and filter differential results
diff_df <- summary(diff_result)
significant_diff <- diff_df #[diff_df$FDR <= fdr, ]

# Step 7: Write outputs
# Write significant differential loops to BEDPE file
message("Writing bedpe differential analysis between ", control_name, " and ", test_name)
output_file <- paste0(output_dir, "/diffloop_", test_name, "over", control_name, "_differential_loops.bedpe")
message("Writing significant differential loops to ", output_file)
write.table(significant_diff, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Save full differential results as RDS for downstream use
rds_file <- paste0(output_dir, "/diffloop_", test_name, "over", control_name, "_differential_results.rds")
message("Saving full differential results to ", rds_file)
saveRDS(diff_result, file = rds_file)

# Print completion message
print("Differential analysis completed.")