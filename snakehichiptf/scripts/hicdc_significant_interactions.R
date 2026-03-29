# Load libraries
library(HiCDCPlus)
library(dplyr)
library(GenomicRanges)
library(optparse)
library(data.table)

# Define command-line arguments
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = NULL, 
              help = "Path to HiC-Pro .allvalidpairs file", metavar = "FILE"),
  make_option(c("-o", "--output_dir"), type = "character", default = "hicdc_output", 
              help = "Output directory [default: %default]", metavar = "DIR"),
  make_option(c("-b", "--bintolen_file"), type = "character", default = NULL, 
              help = "Path to pre-generated bintolen file (e.g., mm10_5000bp_bintolen.txt.gz)", metavar = "FILE"),
  make_option(c("-c", "--chrs"), type = "character", default = NULL, 
              help = "Path to chromosome size file (use first column)", metavar = "FILE"),
  make_option(c("-f", "--fdr"), type = "double", default = 0.01, 
              help = "FDR threshold for significance [default: %default]", metavar = "FDR"),
  make_option(c("-s", "--sample_name"), type = "character", default = NULL, 
              help = "Sample name in the final file  [default: %default]", metavar = "Samplename")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required arguments
if (is.null(opt$input) || is.null(opt$bintolen_file)) {
  print_help(opt_parser)
  stop("Input .hic file and bintolen file are required.", call. = FALSE)
}

# Assign variables
input_file <- opt$input
output_dir <- opt$output_dir
bintolen_file <- opt$bintolen_file
chrom_size_file <- opt$chrs
fdr <- opt$fdr
sample_name <- opt$sample_name

##chrom size
chrom_sizes <- read.table(chrom_size_file, header = FALSE, stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("chr", "length")

get_chrs_custom <- function(chrom_sizes) {
  # chrom_sizes is a data frame with columns 'chr' and 'length'
  return(sort(unique(chrom_sizes$chr)))
}

get_chr_sizes_custom <- function(chrom_sizes, chrs = NULL) {
  if (is.null(chrs)) {
    chrs <- get_chrs_custom(chrom_sizes)
  }
  sizes <- setNames(chrom_sizes$length, chrom_sizes$chr)
  return(sizes[chrs])
}

generate_df_gi_list <- function(df, chrs = NULL, Dthreshold = 2e+06, chrom_sizes) {
  gi_list <- list()
  if (!("chr" %in% colnames(df) & "start" %in% colnames(df))) {
    stop("the df should have columns named 'chr' and 'start'")
  }
  if (!sum(stats::complete.cases(df %>% dplyr::select(.data$chr, .data$start))) == nrow(df)) {
    stop("the df has NULL/NA values for 'chr','start'")
  }
  if (!(sum(is.finite(df$start)) == nrow(df))) {
    stop("the df has NA/NaN/Inf values for 'start'")
  }
  if (is.null(chrs)) {
    chrs <- sort(unique(df$chr))
  }
  chrom_lengths <- get_chr_sizes_custom(chrom_sizes, chrs)
  
  for (chrom in chrs) {
    df_chr <- df %>% dplyr::filter(.data$chr == chrom)
    if (!"end" %in% colnames(df_chr)) {
      df_chr <- df_chr %>% dplyr::mutate(end = c(.data$start[-1], chrom_lengths[chrom]))
    }
    all.regions <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = df_chr$start, end = df_chr$end))
    eff_binsize <- stats::quantile(df_chr$end - df_chr$start, probs = 0.01)
    numbins <- length(all.regions)
    maxbins <- min(round(Dthreshold / eff_binsize), numbins)
    index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
    index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
    gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
    mcols(gi_list[[chrom]])$D <- InteractionSet::pairdist(gi_list[[chrom]])
    gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])$D <= Dthreshold]
  }
  return(gi_list)
}

generate_binned_gi_list <- function(binsize, chrs = NULL, Dthreshold = 2e+06, chrom_sizes) {
  gi_list <- list()
  if (is.null(chrs)) {
    chrs <- get_chrs_custom(chrom_sizes)
  } else {
    chrs <- chrs[chrs %in% get_chrs_custom(chrom_sizes)]
  }
  chrom_lengths <- get_chr_sizes_custom(chrom_sizes, chrs)
  
  for (chrom in chrs) {
    seqlen <- chrom_lengths[chrom]
    numbins <- ceiling(seqlen / binsize)
    maxbins <- min(round(Dthreshold / binsize), numbins)
    all.regions <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = seq(0, (numbins - 1), 1) * binsize, 
                                                                  end = pmin(seq(1, numbins, 1) * binsize, seqlen)))
    index1 <- unlist(lapply(seq(1, numbins, 1), function(x) rep(x, min(maxbins + 1, numbins - x + 1))))
    index2 <- unlist(lapply(seq(1, numbins, 1), function(x) seq(x, min(x + maxbins, numbins), 1)))
    gi_list[[chrom]] <- InteractionSet::GInteractions(index1, index2, all.regions)
    mcols(gi_list[[chrom]])$D <- InteractionSet::pairdist(gi_list[[chrom]])
    gi_list[[chrom]] <- gi_list[[chrom]][mcols(gi_list[[chrom]])$D <= Dthreshold]
  }
  return(gi_list)
}

generate_bintolen_gi_list <- function (bintolen_path, chrs = NULL, Dthreshold = 2e+06, binned = TRUE,binsize = NULL, gen = "Hsapiens", gen_ver = "hg19") {
    input.file.read <- function(filepath) {
        if (grepl("\\.txt.gz$", filepath) | grepl("\\.txt$", 
            filepath)) {
            return(data.table::fread(filepath))
        }
        else if (grepl("\\.rds$", filepath)) {
            return(readRDS(filepath))
        }
        else {
            stop("Can only read in paths ending with .txt,.txt.gz, or .rds")
        }
    }
    if (binned & is.null(binsize) & !is.null(bintolen_path)) {
        bintolen <- input.file.read(bintolen_path)
        bintolen <- bintolen %>% tidyr::separate(.data$bins, 
            c("chr", "start", "end"), "-") %>% dplyr::mutate(start = floor(as.numeric(.data$start)/1000) * 
            1000, end = as.numeric(.data$end))
        binsize <- round(mean(abs(bintolen$end - bintolen$start))/1000) * 
            1000
    }
    else if (binned & is.null(binsize) & is.null(bintolen_path)) {
        stop("If binned, need to specify at least one of binsize and bintolen_path")
    }
    if (!is.null(bintolen_path) & !exists("bintolen")) {
        bintolen <- input.file.read(bintolen_path)
        if (binned) {
            bintolen <- bintolen %>% tidyr::separate(.data$bins, 
                c("chr", "start", "end"), "-") %>% dplyr::mutate(start = floor(as.numeric(.data$start)/1000) * 
                1000, end = as.numeric(.data$end))
        }
        else {
            bintolen <- bintolen %>% tidyr::separate(.data$bins, 
                c("chr", "start", "end"), "-") %>% dplyr::mutate(start = as.numeric(.data$start), 
                end = as.numeric(.data$end))
        }
    }
    if (is.null(chrs)) {
        chrs <- sort(unique(bintolen$chr))
    }
    if (binned) {
        gi_list <- generate_binned_gi_list(binsize, chrs, Dthreshold, chrom_sizes)
    }
    else {
        gi_list <- generate_df_gi_list(bintolen, chrs, Dthreshold, chrom_sizes)
    }
    gi_list <- add_1D_features(gi_list, bintolen, chrs)
    return(gi_list)
}


# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

#generate gi_list instance
gi_list<-generate_bintolen_gi_list(bintolen_path=bintolen_file)

#add .hic counts
gi_list<-add_hicpro_allvalidpairs_counts(gi_list, allvalidpairs_path = input_file)

#expand features for modeling
gi_list<-expand_1D_features(gi_list)

#run HiC-DC+ 
set.seed(1010) #HiC-DC downsamples rows for modeling
gi_list<-HiCDCPlus(gi_list, Dmin = 5000) #HiCDCPlus_parallel runs in parallel across ncores
print("hicdcplus result gi_list")
print(head(gi_list))
save(gi_list, file=paste0(output_dir, paste0("/Hicdcplus.",sample_name,".all_interactions.RData")))

# Export to .hic for visualization
hic_file_out <- paste0(output_dir, "/Hicdcplus.", sample_name, ".all_interactions.hic")

hic_ok <- TRUE
tryCatch({
  hicdc2hic(gi_list, hicfile = hic_file_out, mode = "normcounts")
}, error = function(e) {
  hic_ok <<- FALSE
  message("[WARN] hicdc2hic() failed for sample: ", sample_name)
  message("[WARN] ", conditionMessage(e))
})

# Export to bedpe
# Combine all GInteractions objects and filter for significant interactions
all_gi <- gi_list[[1]]
if (length(gi_list) > 1) {
  for (i in 2:length(gi_list)) {
    all_gi <- c(all_gi, gi_list[[i]])
  }
}

if (!hic_ok) {
  message("[WARN] .hic was not generated, but continuing to export BEDPE/RData outputs.")
}
saveRDS(all_gi, file = paste0(output_dir, paste0("/Hicdcplus.",sample_name,".all_interactions.ginteraction.rds")))

# Convert to data frame, adjust coordinates, and select BEDPE columns
selected_column = c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "D", "counts", "pvalue","qvalue", "mu", "sdev")
bedpe_df <- as.data.frame(all_gi)[, selected_column]
# bedpe_df$start1 <- bedpe_df$start1 - 1
# bedpe_df$start2 <- bedpe_df$start2 - 1
colnames(bedpe_df) = c("chrI", "startI", "endI", "chrJ", "startJ", "endJ", "D", "counts", "pvalue","qvalue", "mu", "sdev")

# Convert to data.table for faster operations
print("exporting bedpe file")
setDT(bedpe_df)
fwrite(bedpe_df, 
       paste0(output_dir, paste0("/Hicdcplus.",sample_name,".all_interactions.bedpe")), sep = "\t", quote = FALSE, 
       row.names = FALSE, col.names = TRUE)

sig_bedpe_df <- bedpe_df[qvalue <= fdr, ]
sig_bedpe_df$start1 <- sig_bedpe_df$start1 - 1
sig_bedpe_df$start2 <- sig_bedpe_df$start2 - 1
fwrite(sig_bedpe_df, 
       paste0(output_dir, paste0("/Hicdcplus.",sample_name,".significant_interactions.bedpe")), sep = "\t", quote = FALSE, 
       row.names = FALSE, col.names = TRUE)

system(paste("gzip -f", paste0(output_dir, paste0("/Hicdcplus.",sample_name,".all_interactions.bedpe"))))
system(paste("gzip -f", paste0(output_dir, paste0("/Hicdcplus.",sample_name,".significant_interactions.bedpe"))))



##Another way using dataframe to 
#interaction_df_list <- lapply(gi_list, function(gi) {
#  as.data.frame(gi)[, c("seqnames1", "start1", "end1", "seqnames2", "start2", "end2", "counts", "qvalue")]
#})
#all_interaction_df <- do.call(rbind, interaction_df_list)
#all_interaction_df$start1 <- all_interaction_df$start1 - 1  # Adjust to 0-based for BEDPE
#all_interaction_df$start2 <- all_interaction_df$start2 - 1
#sig_interaction_df <- all_interaction_df[all_interaction_df$qvalue <= fdr, ]
#write.table(all_interaction_df, paste0(output_dir, "/Hicdcplus.all_interactions.bedpe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(sig_interaction_df, paste0(output_dir, "/Hicdcplus.significant_interactions.bedpe"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)