# Load library
library(HiCDCPlus)
library(Biostrings)
library(optparse)
library(GenomicRanges)
library(dplyr)
library(IRanges)

# Define command-line arguments
option_list <- list(
  make_option(c("-o", "--output_dir"), type = "character", default = "hicdc_features", 
              help = "Output directory for features [default: %default]", metavar = "DIR"),
  make_option(c("-g", "--genome"), type = "character", default = "mm10", 
              help = "Genome fasta [default: %default]", metavar = "GENOME"),
  make_option(c("-r", "--restriction_enzyme"), type = "character", default = "MboI", 
              help = "Restriction enzyme [default: %default]", metavar = "ENZYME"),
  make_option(c("-b", "--bin_size"), type = "integer", default = 5000, 
              help = "Bin size in bp [default: %default]", metavar = "SIZE"),
  make_option(c("-c", "--chrs"), type = "character", default = NULL, 
              help = "Path to chromosome size file (use first column)", metavar = "FILE")
)

# Parse arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check required argument
if (is.null(opt$chrs)) {
  print_help(opt_parser)
  stop("Chromosome size file is required.", call. = FALSE)
}

# Assign variables
output_dir <- opt$output_dir
genome <- opt$genome
restriction_enzyme <- strsplit(opt$restriction_enzyme, split = ",")[[1]]
bin_size <- opt$bin_size
chrom_size_file <- opt$chrs

# Read chromosome sizes (first column)
chrom_sizes <- read.table(chrom_size_file, header = FALSE, stringsAsFactors = FALSE)
chrs <- chrom_sizes$V1  # First column as chromosome names

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

##changed function using customized genome.fa
get_enzyme_cutsites <- function(sig, gen, chrs = NULL) {
  print("Customized get_enzyme_cutsites being used")
  # If chrs not provided, use all chromosomes from gen
  if (is.null(chrs)) {
    chrs <- names(gen)
    chrs <- chrs[!chrs %in% c("chrY", "chrM")]  # Exclude Y and mitochondrial if present
  }
  
  # Expand ambiguous patterns (e.g., "N" to "A", "C", "G", "T")
  while (any(grepl("N", sig))) {
    if (length(sig) <= 1) {
      sig <- unique(c(sub("N", "A", sig), sub("N", "C", sig), 
                      sub("N", "G", sig), sub("N", "T", sig)))
    } else {
      sig_new <- NULL
      for (i in seq_along(sig)) {
        if (grepl("N", sig[i])) {
          sig_new_i <- unique(c(sub("N", "A", sig[i]), sub("N", "C", sig[i]), 
                                sub("N", "G", sig[i]), sub("N", "T", sig[i])))
          sig_new <- c(sig_new, sig_new_i)
        } else {
          sig_new <- unique(c(sig_new, sig[i]))
        }
      }
      sig <- unique(sig_new)
    }
  }
  
  # Find restriction sites for a single chromosome (called from construct_features_chr)
  if (length(chrs) == 1) {
    seq <- gen[[chrs]]
    if (length(sig) <= 1) {
      matches <- matchPattern(sig, seq)
      range_chr <- IRanges(start = start(matches), end = end(matches))
    } else {
      patternobj <- DNAStringSet(sig)
      matches <- matchPDict(patternobj, seq)
      range_chr <- ranges(extractAllMatches(seq, matches))
    }
    if (length(start(range_chr)) <= 0) {
      warning(paste0("Chromosome ", chrs, " does not have any cutsites with pattern set ", 
                     paste(sig, collapse = ","), ". Using full chromosome instead"))
      range_chr <- IRanges(start = 1, end = width(seq))
    }
    enzymeCuts <- sort(GenomicRanges::GRanges(seqnames = chrs, ranges = range_chr))
    return(enzymeCuts)
  }
  
  # For multiple chromosomes (not typically used in this context)
  enzymeCuts <- lapply(chrs, function(x) {
    seq <- gen[[x]]
    if (length(sig) <= 1) {
      matches <- matchPattern(sig, seq)
      range_chr <- IRanges(start = start(matches), end = end(matches))
    } else {
      patternobj <- DNAStringSet(sig)
      matches <- matchPDict(patternobj, seq)
      range_chr <- ranges(extractAllMatches(seq, matches))
    }
    if (length(start(range_chr)) <= 0) {
      warning(paste0("Chromosome ", x, " does not have any cutsites with pattern set ", 
                     paste(sig, collapse = ","), ". Using full chromosome instead"))
      range_chr <- IRanges(start = 1, end = width(seq))
    }
    sort(GenomicRanges::GRanges(seqnames = x, ranges = range_chr))
  })
  enzymeCuts <- suppressWarnings(sort(unlist(as(enzymeCuts, "GRangesList"))))
  return(enzymeCuts)
}

construct_features_chr <- function(chrom, gen, gen_ver = "mm10", sig = "GATC", 
                                   bin_type = "Bins-uniform", binsize = 5000, wg_file = NULL, 
                                   feature_type = "RE-based") {
  print("Customized construct_features_chr being used")
  library(Biostrings)
  library(GenomicRanges)
  library(dplyr)
  library(IRanges)
  
  # Get chromosome size from gen (DNAStringSet)
  if (!chrom %in% names(gen)) stop(paste("Chromosome", chrom, "not found in genome FASTA."))
  genome.chromSize <- width(gen[chrom])
  genome.chromGR <- GenomicRanges::GRanges(seqnames = chrom, 
                                           ranges = IRanges::IRanges(start = 1, end = genome.chromSize))
  
  # Custom getGC using gen for the specific chromosome
  getGC <- function(regions, chrom_seq) {
    # chrom_seq is a single DNAString for the current chromosome
    if (!all(seqnames(regions) == chrom)) {
      stop("All regions must belong to the same chromosome as provided.")
    }
    views <- Views(chrom_seq, ranges(regions))
    gc <- rowSums(alphabetFrequency(views, as.prob = TRUE)[, c("G", "C")])
    return(gc)
  }
  
  # getMap function (unchanged)
  getMap <- function(regions, data) {
    hits <- GenomicRanges::findOverlaps(data, regions, type = "within")
    DT <- data.frame(queryHits = as.data.frame(hits)[[1]], 
                     subjectHits = as.data.frame(hits)[[2]])
    DT <- DT %>% dplyr::mutate(map = data$score[.data$queryHits], 
                               len = GenomicRanges::end(data)[.data$queryHits] - 
                                     GenomicRanges::start(data)[.data$queryHits] + 1) %>% 
                 dplyr::group_by(.data$subjectHits) %>% 
                 dplyr::summarize(avgmap = stats::weighted.mean(map, w = len))
    map <- rep(0, length(regions))
    map[DT$subjectHits] <- DT$avgmap
    return(map)
  }
  
  if (feature_type == "RE-based") {
    msg <- paste0("Using ", chrom, " and cut patterns ", sig)
    message(msg)
    
    enzymeCutsites <- get_enzyme_cutsites(sig, gen, chrom)
    RE_sites <- sort(enzymeCutsites)
    newends <- dplyr::lead(BiocGenerics::start(RE_sites)) - 1
    newends <- c(newends[-length(newends)], genome.chromSize)
    BiocGenerics::end(RE_sites) <- newends
    
    minstart <- min(BiocGenerics::start(RE_sites))
    if (minstart > 1) {
      RE_sites <- GenomeInfoDb::sortSeqlevels(c(GenomicRanges::GRanges(seqnames = chrom, 
                                                                       ranges = IRanges::IRanges(start = 1, end = minstart)), 
                                                 RE_sites))
      RE_sites <- sort(RE_sites)
    }
  } else {
    if (bin_type != "Bins-uniform") 
      stop("feature_type=RE-agnostic requires a fixed binsize.")
    msg <- paste0("Using ", chrom, " RE agnostic")
    message(msg)
  }
  
  if (!is.null(wg_file)) {
    wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
  } else {
    wgdata <- NULL
  }
  
  # Extract chromosome sequence once for efficiency
  chrom_seq <- gen[[chrom]]
  
  if (bin_type == "Bins-RE-sites") {
    endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                   IRanges::IRanges(start = GenomicRanges::start(RE_sites), width = 200))
    BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
    BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
    endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                   IRanges::IRanges(end = GenomicRanges::end(RE_sites), width = 200))
    BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
    BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
    
    gcL <- getGC(regions = endL, chrom_seq = chrom_seq)
    gcR <- getGC(regions = endR, chrom_seq = chrom_seq)
    gc <- (gcL + gcR) / 2
    RE_sites$gc <- gc
    
    if (!is.null(wgdata)) {
      endL <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                     IRanges::IRanges(start = GenomicRanges::start(RE_sites), width = 500))
      BiocGenerics::end(endL) <- pmin(BiocGenerics::end(endL), genome.chromSize)
      BiocGenerics::start(endL) <- pmax(BiocGenerics::start(endL), 1)
      endR <- GenomicRanges::GRanges(seqnames = as.character(GenomicRanges::seqnames(RE_sites)), 
                                     IRanges::IRanges(end = GenomicRanges::end(RE_sites), width = 500))
      BiocGenerics::end(endR) <- pmin(BiocGenerics::end(endR), genome.chromSize)
      BiocGenerics::start(endR) <- pmax(BiocGenerics::start(endR), 1)
      mapL <- getMap(regions = endL, data = wgdata)
      mapR <- getMap(regions = endR, data = wgdata)
      map <- (mapL + mapR) / 2
      RE_sites$map <- map
    } else {
      RE_sites$map <- 0
    }
    
    bintolen <- as.data.frame(RE_sites, stringsAsFactors = FALSE)
    bintolen <- bintolen %>% dplyr::mutate(RE_id = dplyr::row_number())
    if (binsize > 1) {
      bintolen <- bintolen %>% dplyr::mutate(binNumb = floor((dplyr::row_number() - 1) / binsize) + 1)
      bintolen <- bintolen %>% dplyr::group_by(.data$seqnames, .data$binNumb) %>% 
                  dplyr::summarize(start = min(start), end = max(end), gc = mean(gc), map = mean(map), 
                                   RE_id = min(.data$RE_id)) %>% 
                  dplyr::mutate(width = .data$end - .data$start + 1) %>% 
                  dplyr::rename(chr = "seqnames") %>% 
                  dplyr::mutate(bins = paste(.data$chr, .data$start, .data$end, sep = "-"))
    } else {
      bintolen <- bintolen %>% dplyr::rename(chr = "seqnames") %>% 
                  dplyr::mutate(bins = paste(.data$chr, .data$start, .data$end, sep = "-"))
    }
    bintolen <- bintolen %>% dplyr::select(.data$bins, .data$gc, .data$map, .data$width, .data$RE_id)
    bintolen <- as.data.frame(bintolen, stringsAsFactors = FALSE)
    rownames(bintolen) <- bintolen$bins
    return(bintolen)
  }
  
  if (bin_type == "Bins-uniform") {
    bins.chrom <- function(chrom, binsize) {
      cuts_start <- seq(1, genome.chromSize, by = binsize)
      cuts_end <- seq(binsize, genome.chromSize - (genome.chromSize %% binsize) + binsize, by = binsize)
      cuts_end <- c(cuts_end[-length(cuts_end)], genome.chromSize)
      bins <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = cuts_start, end = cuts_end))
      return(bins)
    }
    
    binsGR <- bins.chrom(chrom, binsize)
    names(binsGR) <- paste(as.character(GenomicRanges::seqnames(binsGR)), 
                           GenomicRanges::start(binsGR), GenomicRanges::end(binsGR), sep = "-")
    
    if (!is.null(wg_file)) {
      wgdata <- rtracklayer::import.bw(wg_file, which = genome.chromGR, as = "GRanges")
    } else {
      wgdata <- NULL
    }
    
    if (feature_type == "RE-based") {
      medians <- (GenomicRanges::start(enzymeCutsites) + GenomicRanges::end(enzymeCutsites)) / 2
      FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, 
                                              ranges = IRanges::restrict(IRanges::IRanges(end = medians - 1, width = 500), 
                                                                         start = 1, end = genome.chromSize))
      FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, 
                                              ranges = IRanges::restrict(IRanges::IRanges(start = medians, width = 500), 
                                                                         start = 1, end = genome.chromSize))
      ends <- c(FragmentendsL, FragmentendsR)
      hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
      LR <- data.frame(bins = names(binsGR[hits[[2]]]), 
                       start = GenomicRanges::start(ends[hits[[1]]]), 
                       end = GenomicRanges::end(ends[hits[[1]]]), stringsAsFactors = FALSE)
    } else {
      LR <- data.frame(bins = names(binsGR), start = GenomicRanges::start(binsGR), 
                       end = GenomicRanges::end(binsGR), stringsAsFactors = FALSE)
    }
    
    LRgr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR$start, end = LR$end))
    if (!is.null(wgdata)) {
      LR$map <- getMap(LRgr, wgdata)
    } else {
      LR$map <- 0
    }
    map <- LR %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(map = mean(map))
    
    if (feature_type == "RE-based") {
      ir <- IRanges::IRanges(LR$start, LR$end)
      LR$group <- as.data.frame(GenomicRanges::findOverlaps(ir, IRanges::reduce(ir)))[[2]]
      len <- LR %>% dplyr::group_by(.data$group, .data$bins) %>% 
            dplyr::summarize(start = min(start), end = max(end)) %>% 
            dplyr::mutate(width = .data$end - .data$start + 1) %>% 
            dplyr::group_by(.data$bins) %>% dplyr::summarize(len = sum(.data$width))
    } else {
      len <- NULL
    }
    
    if (feature_type == "RE-based") {
      FragmentendsL <- GenomicRanges::GRanges(seqnames = chrom, 
                                              ranges = IRanges::IRanges(end = medians - 1, width = 200))
      FragmentendsR <- GenomicRanges::GRanges(seqnames = chrom, 
                                              ranges = IRanges::IRanges(start = medians, width = 200))
      ends <- c(FragmentendsL, FragmentendsR)
      hits <- as.data.frame(GenomicRanges::findOverlaps(ends, binsGR, type = "within", select = "all"))
      LR2 <- data.frame(bins = names(binsGR[hits[[2]]]), 
                        start = GenomicRanges::start(ends[hits[[1]]]), 
                        end = GenomicRanges::end(ends[hits[[1]]]), stringsAsFactors = FALSE)
    } else {
      LR2 <- LR %>% dplyr::select(.data$bins, .data$start, .data$end)
    }
    
    LR2gr <- GenomicRanges::GRanges(chrom, IRanges::IRanges(start = LR2$start, end = LR2$end))
    LR2$gc <- getGC(LR2gr, chrom_seq)
    gc <- LR2 %>% dplyr::group_by(.data$bins) %>% dplyr::summarize(gc = mean(gc))
    
    if (!is.null(len)) {
      LR <- dplyr::left_join(dplyr::left_join(gc, map), len)
    } else {
      LR <- dplyr::left_join(gc, map)
    }
    
    bintolen <- as.data.frame(LR)
    if (sum(bintolen$map) == 0) {
      bintolen <- bintolen %>% dplyr::select(-.data$map)
    }
    rownames(bintolen) <- bintolen$bins
    
    allbins <- data.frame(bins = names(binsGR), stringsAsFactors = FALSE)
    bintolen <- suppressWarnings(dplyr::left_join(allbins, bintolen) %>% 
                                 tidyr::replace_na(list(gc = 0, map = 0, len = 0)))
    return(bintolen)
  }
}

construct_features <- function(output_path, gen_fa, gen_ver = "mm10", sig = "GATC", 
                               bin_type = "Bins-uniform", binsize = 5000, wg_file = NULL, 
                               chrs = NULL, feature_type = "RE-based") {
  print("Customized construct_features being used")
  gen <- readDNAStringSet(gen_fa)
  
  # If chrs not provided, use all chromosomes from the FASTA
  if (is.null(chrs)) {
    chrs <- names(gen)
  }
  
  # Process each chromosome
  bintolen <- lapply(chrs, function(x) data.frame(stringsAsFactors = FALSE))
  names(bintolen) <- chrs
  for (chrom in names(bintolen)) {
    bintolen[[chrom]] <- construct_features_chr(
      chrom = chrom,
      gen = gen,  # Pass the DNAStringSet object
      gen_ver = gen_ver,
      sig = sig,
      bin_type = bin_type,
      binsize = binsize,
      wg_file = wg_file,
      feature_type = feature_type
    )
  }
  
  # Combine and write output
  bintolen_combined <- suppressWarnings(dplyr::bind_rows(bintolen))
  bintolenoutput <- path.expand(paste0(output_path, ".bintolen.txt.gz"))
  bintolenoutputdir <- dirname(bintolenoutput)
  if (!dir.exists(bintolenoutputdir)) {
    dir.create(bintolenoutputdir, showWarnings = FALSE, recursive = TRUE, mode = "0777")
  }
  data.table::fwrite(bintolen_combined, bintolenoutput, row.names = FALSE, quote = FALSE, sep = "\t")
  
  return(bintolenoutput)
}

print(paste0("restriction_enzyme: ",restriction_enzyme))
# Generate genomic features
features_path <- paste0(output_dir, "/", "restriction_enzyme.hicdcplus")
construct_features(
  output_path = features_path,
  gen = genome,
  gen_ver = "",
  sig = restriction_enzyme,
  bin_type = "Bins-uniform",
  binsize = bin_size,
  chrs = chrs
)

cat("Genomic features generated at:", features_path, "\n")
cat("Bintolen file:", paste0(features_path, ".bintolen.txt.gz"), "\n")
