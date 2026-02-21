options("scipen"=999)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

INDIR = args[1]
SET = args[2]
FDR = if(-log10(as.numeric(args[3])) %% 1 == 0) as.integer(-log10(as.numeric(args[3]))) else -log10(as.numeric(args[3]))
RESOLUTION = as.numeric(args[4])

# Input file for significant peaks
inf = paste(INDIR, SET, '.', as.character(FDR), '.peaks', sep='')
peaks_raw = fread(inf, header=TRUE, stringsAsFactors=FALSE)

# Define the final column names
final_colnames <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', 'count', 
                    'expected', 'fdr', 'ClusterLabel', 'ClusterSize', 'ClusterType', 
                    'ClusterNegLog10P', 'ClusterSummit')

# Process significant interactions
peaks = as.data.table(subset(peaks_raw, ClusterType != 'Singleton' | 
                            (ClusterType == 'Singleton' & fdr < 1e-4)))
if (nrow(peaks) > 0) {
    peaks[, summit := 1 * (fdr == min(fdr)), by = lab]
    peaks$summit[peaks$ClusterType == 'Singleton'] = 1

    singleton_labs = peaks$lab[peaks$ClusterType == 'Singleton']
    peaks$lab[peaks$ClusterType == 'Singleton'] = paste(singleton_labs, 1:length(singleton_labs), sep='')

    peaks$bin1_end = peaks$bin1_mid + RESOLUTION
    peaks$bin2_end = peaks$bin2_mid + RESOLUTION
    peaks_final = subset(peaks, select = c("chr", "bin1_mid", "bin1_end", "chr", 
                                          "bin2_mid", "bin2_end", "count", "expected2", 
                                          "fdr", "lab", "ClusterSize", "ClusterType", 
                                          "NegLog10P", "summit"))
    colnames(peaks_final) = final_colnames
} else {
    peaks_final <- data.frame(
        chr1 = character(0),
        start1 = integer(0),
        end1 = integer(0),
        chr2 = character(0),
        start2 = integer(0),
        end2 = integer(0),
        count = integer(0),
        expected = numeric(0),
        fdr = numeric(0),
        ClusterLabel = character(0),
        ClusterSize = integer(0),
        ClusterType = character(0),
        ClusterNegLog10P = numeric(0),
        ClusterSummit = integer(0)
    )
}

# Output significant interactions with # header
fout_sig = paste(INDIR, "Maps.", SET, '.sig3Dinteractions.bedpe', sep='')
print(paste("Significant interactions filename: ", fout_sig))
write(paste("#", paste(final_colnames, collapse="\t"), sep=""), file=fout_sig)
fwrite(peaks_final, fout_sig, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t', append=TRUE)

# Input file for all peaks
all_inf = paste(INDIR, SET, '.all.peaks', sep='')

if (file.exists(all_inf)) {
    all_peaks_raw = fread(all_inf, header=TRUE, stringsAsFactors=FALSE)

    # Process all interactions
    all_peaks = as.data.table(all_peaks_raw)
    if (nrow(all_peaks) > 0) {
        all_peaks[, summit := 1 * (fdr == min(fdr)), by = lab]
        all_peaks$bin1_end = all_peaks$bin1_mid + RESOLUTION
        all_peaks$bin2_end = all_peaks$bin2_mid + RESOLUTION
        all_peaks_final = subset(all_peaks, select = c("chr", "bin1_mid", "bin1_end", "chr", 
                                                    "bin2_mid", "bin2_end", "count", "expected2", 
                                                    "fdr", "lab", "ClusterSize", "ClusterType", 
                                                    "NegLog10P", "summit"))
        colnames(all_peaks_final) = final_colnames
    } else {
        all_peaks_final <- data.frame(
            chr1 = character(0),
            start1 = integer(0),
            end1 = integer(0),
            chr2 = character(0),
            start2 = integer(0),
            end2 = integer(0),
            count = integer(0),
            expected = numeric(0),
            fdr = numeric(0),
            ClusterLabel = character(0),
            ClusterSize = integer(0),
            ClusterType = character(0),
            ClusterNegLog10P = numeric(0),
            ClusterSummit = integer(0)
        )
    }

    # Output all interactions with # header
    fout_all = paste(INDIR, "Maps.", SET, '.all3Dinteractions.bedpe', sep='')
    print(paste("All interactions filename: ", fout_all))
    write(paste("#", paste(final_colnames, collapse="\t"), sep=""), file=fout_all)
    fwrite(all_peaks_final, fout_all, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t', append=TRUE)
}