library(ggplot2)
  
# Function to create quantile dot plot
create_quantile_dot_plot <- function(file_pattern, output_file, ref_sample_order = NULL, sample_list = NULL, plot_title = "") {
  # Load required library
  # Get list of quantile files
  quantile_files <- Sys.glob(file_pattern)
  
  # Check if files are found
  if (length(quantile_files) == 0) {
    stop("No files found matching the pattern")
  }
  
  # Initialize data frame to store all quantiles
  all_data <- data.frame(Quantile = character(), Value = numeric(), Sample = character(), stringsAsFactors = FALSE)
  
  # Read quantile files
  for (file in quantile_files) {
    if (file.exists(file)) {
      quantiles <- read.csv(file, stringsAsFactors = FALSE)
      # Extract sample name from directory
      sample <- basename(dirname(dirname(file)))
      if(length(grep("atac",file_pattern)) !=0){sample <- basename(dirname(file))}
      quantiles$Sample <- sample
      all_data <- rbind(all_data, quantiles)
    } else {
      warning(paste("File not found:", file))
    }
  }
  
  # Validate data
  if (nrow(all_data) == 0) {
    stop("No quantile data found in any sample directories")
  }
  
  # Filter by sample_list if provided (for HiCHIP)
  #if (!is.null(sample_list)) {
  #  missing_samples <- setdiff(sample_list, all_data$Sample)
  #  if (length(missing_samples) > 0) {
  #    warning(paste("The following samples were not found in the data:", 
  #                  paste(missing_samples, collapse = ", ")))
  #  }
  #  all_data <- all_data[all_data$Sample %in% sample_list, ]
  #  if (nrow(all_data) == 0) {
  #    stop("No data remains after filtering by sample list")
  #  }
  #}
  
  # Extract second item from sample names and handle duplicates
  sample_names <- unique(all_data$Sample)
  second_items <- sapply(strsplit(sample_names, "_"), function(x) if (length(x) >= 2) x[2] else x[1])
  # Create unique labels for duplicates
  unique_labels <- character(length(sample_names))
  name_counts <- table(second_items)
  for (i in seq_along(second_items)) {
    name <- second_items[i]
    if (name_counts[name] > 1) {
      indices <- which(second_items == name)
      unique_labels[i] <- paste0(name, "_", which(indices == i))
    } else {
      unique_labels[i] <- name
    }
  }
  # Map original sample names to new labels
  label_map <- setNames(unique_labels, sample_names)
  all_data$Sample_Label <- label_map[all_data$Sample]
  
  # Check for unmatched samples in ATAC when ref_sample_order is provided
  if (!is.null(ref_sample_order)) {
    unmatched_samples <- setdiff(unique_labels, ref_sample_order)
    if (length(unmatched_samples) > 0) {
      warning(paste("The following samples in the dataset are not in the reference sample order:", 
                    paste(unmatched_samples, collapse = ", ")))
    }
  }
  
  # Filter out 50th quantile
  all_data <- subset(all_data, Quantile != "25th")
  all_data <- subset(all_data, Quantile != "50th")
    
  # Ensure Quantile is ordered correctly
  all_data$Quantile <- factor(all_data$Quantile, levels = c("75th", "90th", "95th", "99th"))
  print(levels(all_data$Quantile))

  # Order samples
  if (!is.null(ref_sample_order)) {
    print("Using provided reference sample order.")
    # Use provided reference sample order, keeping only samples present in current data
    valid_ref_order <- ref_sample_order[ref_sample_order %in% unique_labels]
    all_data$Sample_Label <- factor(all_data$Sample_Label, levels = valid_ref_order)
    print(paste("Using reference sample order with", length(valid_ref_order), "samples."))
    print(levels(all_data$Sample_Label))
  } else {
    print("No reference sample order provided; ordering by 99th quantile values.")
    # Order by 99th quantile value (ascending)
    quantile_99 <- subset(all_data, Quantile == "99th")
    sample_order <- quantile_99$Sample[order(quantile_99$Value, decreasing = FALSE)]
    print(all_data)
    all_data$Sample_Label <- factor(all_data$Sample_Label, levels = label_map[sample_order])
  }
  
  # Assign shapes cyclically (3 shapes: circle=16, triangle=17, square=15)
  shapes <- c(16, 17, 15)  # Circle, triangle, square
  all_data$Shape <- shapes[(as.numeric(all_data$Sample_Label) - 1) %% 3 + 1]
  
  # Create position for side-by-side dots within each quantile
  all_data$Position <- as.numeric(all_data$Sample_Label)  # Numeric position based on sample order
  
  # Create plot
  p <- ggplot(all_data, aes(x = Position, y = Log1p_Value, color = Sample_Label, shape = Sample_Label)) +
    geom_point(size = 2) +
    facet_grid(. ~ Quantile, scales = "free_x", space = "free_x") +
    scale_x_continuous(
      breaks = seq_along(levels(all_data$Sample_Label)),
      labels = NULL  # Remove x-axis labels
    ) +
    scale_shape_manual(values = shapes[rep(1:3, length.out = length(levels(all_data$Sample_Label)))]) +
    theme_classic() +
    labs(
      title = plot_title,
      x = NULL,  # No x-axis title
      y = "log1p(Counts)",
      color = "Samples",
      shape = "Samples"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_blank(),  # Remove x-axis tick labels
      axis.ticks.x = element_blank(),  # Remove x-axis ticks
      panel.grid.major = element_blank(),  # Remove major grid
      panel.grid.minor = element_blank(),  # Remove minor grid
      panel.border = element_blank(),  # Remove panel borders
      plot.background = element_rect(color = "black", fill = NA, size = 0.5),  # Frame around entire plot
      legend.position = c(0.01, 0.99),  # Legend inside plot, top-left
      legend.justification = c(0, 1),   # Align to top-left
      legend.text = element_text(size = 5, family = "Helvetica"),
      legend.background = element_rect(fill = NA, color = NA),
      strip.text = element_text(size = 12, family = "Helvetica"),
      axis.text.y = element_text(size = 12, family = "Helvetica", colour = "black")
    )
  
  # Save plot
  ggsave(output_file, plot = p, width = 12, height = 6, dpi = 300)
  
  # Print summary for verification
  print(summary(all_data))
  
  # Check quantile similarity (relative difference < 10% on log1p values)
  rel_diff <- aggregate(Log1p_Value ~ Quantile, data = all_data, function(x) (max(x) - min(x)) / mean(x))
  if (any(rel_diff$Log1p_Value > 0.1)) {
    warning("Log1p quantile differences exceed 10% for some quantiles")
  }
  print(rel_diff)
  
  # Return sample order for use in subsequent plots
  return(levels(all_data$Sample_Label))
}

# Define the list of HiCHIP samples to include
hichip_samples <- c(
  "Human_acinar26.6_unknown_WT_unknown_HiCHIP_standard",
  "Human_alpha24.4_unknown_WT_unknown_HiCHIP_standard",
  "Human_Astrocyte_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_Bcellprecursorleukemia_unknown_WT_unknown_HiCHIP_standard",
  "Human_beta24.4_unknown_WT_unknown_HiCHIP_standard",
  "Human_BladderEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_BronchialTrachealEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_ColonicEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_delta24.4_unknown_WT_unknown_HiCHIP_standard",
  "Human_duct35.8_unknown_WT_unknown_HiCHIP_standard",
  "Human_EndometrialEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_GM12878_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_GP5d_unknown_WT_unknown_HiCHIP_standard",
  "Human_Jurkat_unknown_WT_unknown_HiCHIP_standard",
  "Human_MammaryEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_NaiveT_unknown_WT_unknown_HiCHIP_standard",
  "Human_OCILy7_unknown_WT_unknown_HiCHIP_standard",
  "Human_ProstateEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  "Human_ProximalTubularEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR",
  #"Human_Caudate_unknow_WT_unknow_HiCHIP_standard_GSM4441834",
  #"Human_Caudate_unknow_WT_unknow_HiCHIP_standard_GSM4441835",
  #"Human_Hippocampus_unknow_WT_unknow_HiCHIP_standard_GSM4441836",
  #"Human_Hippocampus_unknow_WT_unknow_HiCHIP_standard_GSM4441837",
  #"Human_MiddleFrontalGyrus_unknow_WT_unknow_HiCHIP_standard_GSM4441830",
  #"Human_MiddleFrontalGyrus_unknow_WT_unknow_HiCHIP_standard_GSM4441833",
  #"Human_ParietalLobe_unknow_WT_unknow_HiCHIP_standard_GSM4441838",
  #"Human_ParietalLobe_unknow_WT_unknow_HiCHIP_standard_GSM4441839",
  #"Human_SubstantiaNigra_unknow_WT_unknow_HiCHIP_standard_GSM4441831",
  #"Human_SubstantiaNigra_unknow_WT_unknow_HiCHIP_standard_GSM4441832",
  #"Human_SuperiorandMiddleTemporalGyri_unknow_WT_unknow_HiCHIP_standard_GSM4441840",
  #"Human_SuperiorandMiddleTemporalGyri_unknow_WT_unknow_HiCHIP_standard_GSM4441841",
  "Human_Th17_unknown_WT_unknown_HiCHIP_standard",
  "Human_ThyroidEpithelial_unknown_WT_unknown_HiCHIP_standard_mergedSRR"
)
hichip_samples=NULL

# Get command line arguments
if(FALSE){
if (length(grep("hichip", args[1], ignore.case = TRUE)) == 10) {
  # HiCHIP data with filtered (>0) samples
  hichip_file_pattern <- args[2] #"/storage2/fs1/syidan/Active/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/*/hichip_matrix/quantiles_nonzero.csv"
  hichip_sample_order <- create_quantile_dot_plot(
    file_pattern = hichip_file_pattern,
    output_file = paste0(file.path(args[1],"quantile_dot_plot_hichip.pdf")),
    sample_list = hichip_samples,
    plot_title = "HiCHIP Quantile Plot"
  )
  print(hichip_sample_order)

  # HiCHIP data with filtered (>0) samples
  hichip_file_pattern <- args[3] #"/storage2/fs1/syidan/Active/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/*/hichip_matrix/quantiles_all.csv"
  create_quantile_dot_plot(
    file_pattern = hichip_file_pattern,
    output_file = paste0(file.path(args[1],"quantile_dot_plot_hichip_over1.pdf")),
    sample_list = hichip_samples,
    plot_title = "HiCHIP Quantile Plot"
  )
  create_quantile_dot_plot(
    file_pattern = hichip_file_pattern,
    output_file = paste0(file.path(args[1],"quantile_dot_plot_hichip_over1_order.pdf")),
    sample_list = hichip_samples,
    ref_sample_order = hichip_sample_order,
    plot_title = "HiCHIP Quantile Plot"
  )
} else if (length(grep("atac", args[1], ignore.case = TRUE)) == 10) {
  # ATAC data filtered (>0), using HiCHIP sample order
  atac_file_pattern <- args[2] #"/storage2/fs1/syidan/Active/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/*/genomic_features/atac_quantiles_nonzero.csv"
  create_quantile_dot_plot(
    file_pattern = atac_file_pattern,
    output_file = paste0(file.path(args[1],"quantile_dot_plot_atac.pdf")),
    #ref_sample_order = hichip_sample_order,
    plot_title = "ATAC Quantile Plot"
  )
  }
}

# HiCHIP data with filtered (>0) samples
args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1] #"/storage2/fs1/syidan/Active/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/quantiles"
hichip_file_pattern <- args[2] #"/storage2/fs1/syidan/Active/Data/Processed/HiCHIP_GSE_softlink/HiCHIP/merge_data/*/hichip_matrix/quantiles_nonzero.csv"
filename <- args[3] #"hichip"
#print(output_dir)
#print(hichip_file_pattern)
hichip_sample_order <- create_quantile_dot_plot(
  file_pattern = hichip_file_pattern,
  output_file = paste0(output_dir,"/quantile_dot_plot_",filename,".pdf"),
  sample_list = hichip_samples,
  plot_title = paste0(filename, "Quantile Plot")
)
#print(hichip_sample_order)

#write("Plotting completed successfully.\n", file=paste0(file.path(output_dir,"quantile_complete.txt")))