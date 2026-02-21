#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(HiCDCPlus)
})

# ----------------------------
# Command line arguments
# ----------------------------
option_list <- list(
  make_option("--samplesheet", type="character", help="TSV with columns: path, condition"),
  make_option("--output_dir", type="character", help="Output directory"),
  make_option("--fdr", type="double", default=0.05, help="FDR cutoff for union index (default=0.05)"),
  make_option("--binsize", type="integer", default=5000),
  make_option("--fitType", type="character", default="mean")
)

opt <- parse_args(OptionParser(option_list=option_list))

if (is.null(opt$samplesheet) || is.null(opt$output_dir)) {
  stop("Must provide --samplesheet and --output_dir")
}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive=TRUE)

if (!grepl("/$", opt$output_dir)) {opt$output_dir <- paste0(opt$output_dir, "/")}


# ----------------------------
# Read samplesheet
# ----------------------------
ss <- read.delim(opt$samplesheet, header=TRUE, stringsAsFactors=FALSE)
ss$path <- file.path(
  opt$output_dir, "samples_link",
  sub("\\.bedpe\\.gz$", ".txt.gz",
      basename(ss$samplename))
)
filter_file <- file.path(opt$output_dir, paste(paste(unique(ss$condition), collapse = "_"), "union_analysis_indices.txt.gz",sep="_"))

if (!all(c("samplename","condition") %in% names(ss))) {
  stop("Samplesheet must contain columns: samplename, condition")
}

# Check files exist
if (any(!file.exists(ss$samplename))) {
  stop("Some files in samplesheet do not exist.")
}

# Link file paths to output directory
dir.create(file.path(opt$output_dir, "samples_link"), showWarnings = FALSE)
file.symlink(
  from = ss[,1],
  to   = file.path(
           opt$output_dir, "samples_link",
           sub("\\.bedpe\\.gz$", ".txt.gz",
               basename(ss[,1]))
         )
)

# ----------------------------
# Step 1: Build union index
# ----------------------------
message("Building union index...")

extract_index <- function(f, fdr_cutoff) { 
  dt <- fread(f, showProgress=FALSE)
  dt[qvalue <= fdr_cutoff, .(chr=chrI, startI=startI, startJ=startJ)] 
  } 
index_list <- lapply(ss$path, extract_index, fdr_cutoff=opt$fdr)
indexfile <- unique(rbindlist(index_list))
setorder(indexfile, chr, startI, startJ)

fwrite(indexfile, filter_file, sep="\t", quote=FALSE)

message("Union index size: ", nrow(indexfile))

# ----------------------------
# Step 2: Construct input_paths
# ----------------------------
input_paths <- split(ss$path, ss$condition)

# Ensure each condition has >=2 replicates
rep_counts <- table(ss$condition)
if (any(rep_counts < 2)) {
  stop("Each condition must have at least 2 replicates for hicdcdiff.")
}

# ----------------------------
# Step 3: Run hicdcdiff
# ----------------------------
message("Running hicdcdiff...")

hicdcdiff(
  input_paths = input_paths,
  filter_file = filter_file,
  output_path = opt$output_dir,
  fitType = opt$fitType,
  binsize = opt$binsize,
  diagnostics = TRUE
)

writeLines(capture.output(sessionInfo()), file.path(opt$output_dir, "hicdcdiff_sessionInfo.txt"))
message("Differential analysis completed.")





##Not used
##Test code
# opt=list()
# opt$samplesheet = read.delim("/home/syidan/syidan/Projects/SnakeHichip/snakehichip/data/example_data/samplesheet.tsv", header=TRUE, stringsAsFactors=FALSE)
# opt$output_dir = "~/syidan/Projects/SnakeHichipResult/ProcessedData/single_enhancer/"
# opt$fdr = 0.05
# opt$binsize = 5000
# opt$fitType = "mean"


# idx <- rbindlist(lapply(ss$samplename, function(p){
#   dt <- fread(p)
#   # dt[, c("startI", "startJ") := .(startI, startJ)]
#   dt[, c("startI", "startJ") := .(startI + 1L, startJ + 1L)]
#   dst <- file.path(opt$output_dir, "samples_link", sub("bedpe\\.gz$", "txt.gz", basename(p)))
#   fwrite(dt, dst, sep = "\t", row.names=FALSE, quote=FALSE, compress = "gzip")
#   dt <- dt[qvalue <= 0.05, .(chr = chrI, startI = startI, startJ = startJ)]
#   unique(dt)
# }))
# indexfile <- unique(idx)
