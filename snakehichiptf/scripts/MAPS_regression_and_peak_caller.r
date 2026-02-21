#!/usr/bin/env Rscript

library(VGAM)
library(MASS)
options(warn=-1)

### Constants
chroms = NULL
runs = c(1)
RESOLUTION = NULL
COUNT_CUTOFF = 5  # Default from original
RATIO_CUTOFF = 2.0
GAP = 15000
FDR = c(2)

### Argument Parsing
args <- commandArgs(trailingOnly=TRUE)
fltr = data.frame(chr='chrNONE', bin=-1)

if (length(args) < 6 || length(args) > 10) {
    print('Wrong number of arguments. Stopping.')
    print('Arguments needed (in this order): INFDIR, SET, RESOLUTION, chroms_file, sex_chroms, regression_type, [COUNT_CUTOFF, RATIO_CUTOFF, FDR, FILTER].')
    print('Optional: COUNT_CUTOFF, RATIO_CUTOFF, FDR, FILTER (use "None" for no filter).')
    print(paste('Number of arguments entered:', length(args)))
    print('Arguments entered:')
    print(args)
    quit()
} else {
    print(args)
    INFDIR = args[1]
    SET = args[2]
    RESOLUTION = as.integer(args[3])
    n_chroms = as.integer(readLines(args[4], n=1))
    sex_chroms = args[5]
    if (grepl('XY', sex_chroms)) {
        chroms = paste('chr', seq(1, n_chroms, 1), sep='')
        chroms = c(chroms, 'chrX', 'chrY')
    } else if (grepl('X', sex_chroms)) {
        chroms = paste('chr', seq(1, n_chroms, 1), sep='')
        chroms = c(chroms, 'chrX')
    } else if (grepl('Y', sex_chroms)) {
        chroms = paste('chr', seq(1, n_chroms, 1), sep='')
        chroms = c(chroms, 'chrY')
    } else {
        chroms = paste('chr', seq(1, n_chroms, 1), sep='')
    }
    REG_TYPE = args[6]
    if (!(REG_TYPE %in% c('pospoisson', 'negbinom'))) {
        print(paste('Wrong regression choice. Your choice:', REG_TYPE, '. Available choices: pospoisson or negbinom'))
        quit()
    }
    if (length(args) >= 7) COUNT_CUTOFF = as.numeric(args[7])
    if (length(args) >= 8) RATIO_CUTOFF = as.numeric(args[8])
    if (length(args) >= 9) FDR = -log10(as.numeric(strsplit(args[9], ",")[[1]]))
    if (length(args) == 10 && args[10] != 'None') {
        fltr = read.table(args[10], header=TRUE)
    }
}

print('Filter used (if any):')
if (length(args) < 10 || args[10] == 'None') print('None') else print(fltr)

### Loading Data
mm_combined_and = data.frame()
mm_combined_xor = data.frame()
outf_names = c()
for (i in chroms) {
    for (j in c('.and', '.xor')) {
        print(paste('Loading chromosome ', i, ' ', j, sep=''))
        inf_name = paste(INFDIR, 'reg_raw.', i, '.', SET, sep='')
        file_path = paste(inf_name, j, sep='')
        outf_names = c(outf_names, paste(inf_name, j, '.MAPS2_', REG_TYPE, sep=''))
        if (file.exists(file_path)) {
            mm = read.table(file_path, header=TRUE)
            if (nrow(mm) > 0) {
                mm$chr = i  # Add chr column as in original
                mm = subset(mm, dist > 1)
                mm = subset(mm, !(mm$chr %in% fltr$chr & (mm$bin1_mid %in% fltr$bin | mm$bin2_mid %in% fltr$bin)))
                if (nrow(mm) == 0) print(paste("Warning: No data remains after filtering", file_path))
            } else {
                print(paste("Warning: ", file_path, " is empty", sep=''))
                mm = data.frame()
            }
        } else {
            print(paste("Warning: ", file_path, " does not exist", sep=''))
            mm = data.frame()
        }
        if (j == '.and') {
            mm_combined_and = rbind(mm_combined_and, mm)
        } else if (j == '.xor') {
            mm_combined_xor = rbind(mm_combined_xor, mm)
        }
    }
}

dataset_length = nrow(mm_combined_and) + nrow(mm_combined_xor)
if (dataset_length == 0) {
    print("Error: No valid data loaded from .and or .xor files")
    quit()
}

### Regression Functions
pospoisson_regression <- function(mm, dataset_length) {
    if (nrow(mm) == 0) {
        print("Warning: Empty data frame, skipping regression")
        return(mm)
    }
    fit <- tryCatch(
        vglm(count ~ logl + loggc + logm + logdist + logShortCount, family=pospoisson(), data=mm),
        error = function(e) {
            print(paste("Warning: First regression failed -", e$message, "- using default expected values"))
            return(NULL)
        }
    )
    if (is.null(fit)) mm$expected <- 1 else mm$expected <- fitted(fit)
    mm$p_val = ppois(mm$count, mm$expected, lower.tail=FALSE, log.p=FALSE) / ppois(0, mm$expected, lower.tail=FALSE, log.p=FALSE)
    m1 = mm[mm$p_val > 1/length(mm$p_val), ]
    fit2 <- tryCatch(
        vglm(count ~ logl + loggc + logm + logdist + logShortCount, family=pospoisson(), data=m1),
        error = function(e) {
            print(paste("Warning: Second regression failed -", e$message, "- using first fit instead"))
            return(NULL)
        }
    )
    if (is.null(fit2)) {
        mm$expected2 <- mm$expected
    } else {
        coeff <- round(coef(fit2), 10)
        mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm +
                                  coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10)
        mm$expected2 <- mm$expected2 / (1 - exp(-mm$expected2))
    }
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = ppois(mm$count, mm$expected2, lower.tail=FALSE, log.p=FALSE) / ppois(0, mm$expected2, lower.tail=FALSE, log.p=FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * dataset_length
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

negbinom_regression <- function(mm, dataset_length) {
    if (nrow(mm) == 0) {
        print("Warning: Empty data frame, skipping regression")
        return(mm)
    }
    fit <- tryCatch(
        glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data=mm),
        error = function(e) {
            print(paste("Warning: First regression failed -", e$message, "- using default expected values"))
            return(NULL)
        }
    )
    if (is.null(fit)) mm$expected <- 1 else mm$expected <- fitted(fit)
    sze = ifelse(is.null(fit), 1, fit$theta)
    mm$p_val = pnbinom(mm$count, mu=mm$expected, size=sze, lower.tail=FALSE)
    m1 = mm[mm$p_val > (1 / length(mm$p_val)), ]
    fit2 <- tryCatch(
        glm.nb(count ~ logl + loggc + logm + logdist + logShortCount, data=m1),
        error = function(e) {
            print(paste("Warning: Second regression failed -", e$message, "- using first fit instead"))
            return(NULL)
        }
    )
    if (is.null(fit2)) {
        mm$expected2 <- mm$expected
    } else {
        coeff <- round(fit2$coefficients, 10)
        mm$expected2 <- round(exp(coeff[1] + coeff[2]*mm$logl + coeff[3]*mm$loggc + coeff[4]*mm$logm +
                                  coeff[5]*mm$logdist + coeff[6]*mm$logShortCount), 10)
    }
    sze2 = ifelse(is.null(fit2), sze, fit2$theta)
    mm$ratio2 <- mm$count / mm$expected2
    mm$p_val_reg2 = pnbinom(mm$count, mu=mm$expected2, size=sze2, lower.tail=FALSE)
    mm$p_bonferroni = mm$p_val_reg2 * dataset_length
    mm$fdr <- p.adjust(mm$p_val_reg2, method='fdr')
    return(mm)
}

### Summary and Peak Functions (unchanged)
do_summaries <- function(peaks_and, peaks_xor, peaks, fraction, r) {
    if (ncol(peaks_and) == 0) peaks_and = data.frame(count=NA, dist=NA, p_val_reg2=NA, fdr=NA)
    if (ncol(peaks_xor) == 0) peaks_xor = data.frame(count=NA, dist=NA, p_val_reg2=NA, fdr=NA)
    if (ncol(peaks_and) == 0 & ncol(peaks_xor) == 0) peaks = data.frame(count=NA, dist=NA, p_val_reg2=NA, fdr=NA)
    summary_one_fdr_val = data.frame(
        run = r, log10_fdr_cutoff = fdr_cutoff, singleton_fraction = fraction,
        AND_size = length(peaks_and$count), AND_mean_dist = mean(peaks_and$dist, na.rm=TRUE) * RESOLUTION,
        AND_median_dist = median(peaks_and$dist, na.rm=TRUE) * RESOLUTION,
        AND_min_count = min(peaks_and$count, na.rm=TRUE), AND_max_pval = max(peaks_and$p_val_reg2, na.rm=TRUE),
        AND_max_fdr = max(peaks_and$fdr, na.rm=TRUE),
        XOR_size = length(peaks_xor$count), XOR_mean_dist = mean(peaks_xor$dist, na.rm=TRUE) * RESOLUTION,
        XOR_median_dist = median(peaks_xor$dist, na.rm=TRUE) * RESOLUTION,
        XOR_min_count = min(peaks_xor$count, na.rm=TRUE), XOR_max_pval = max(peaks_xor$p_val_reg2, na.rm=TRUE),
        XOR_max_fdr = max(peaks_xor$fdr, na.rm=TRUE),
        size = length(peaks$count), mean_dist = mean(peaks$dist, na.rm=TRUE) * RESOLUTION,
        median_dist = median(peaks$dist, na.rm=TRUE) * RESOLUTION,
        min_count = min(peaks$count, na.rm=TRUE), max_pval = max(peaks$p_val_reg2, na.rm=TRUE),
        max_fdr = max(peaks$fdr, na.rm=TRUE)
    )
    return(summary_one_fdr_val)
}

label_peaks <- function(df) {
    if (nrow(df) == 0) {
        print("No peaks to label.")
        return(df)
    }
    chroms = unique(df$chr)
    print('Chromosomes with potential interactions:')
    print(chroms)
    final = data.frame()
    for (CHR in chroms) {
        y = df[df$chr == CHR, ]
        y$p_val_reg2[y$p_val_reg2 == 0] = 1111111
        y$p_val_reg2[y$p_val_reg2 == 1111111] = min(y$p_val_reg2)
        for (i in 1:nrow(y)) {
            z <- y[abs(y$bin1_mid - y$bin1_mid[i]) <= GAP & abs(y$bin2_mid - y$bin2_mid[i]) <= GAP, ]
            y$CountNei[i] <- nrow(z)
        }
        u <- y[y$CountNei == 1, ]
        v <- y[y$CountNei >= 2, ]
        out <- NULL
        if (nrow(u) > 0) {
            u$label <- 0
            u$NegLog10P <- -log10(u$p_val_reg2)
            u$ClusterSize <- 1
            out <- rbind(out, u)
        }
        if (nrow(v) > 0) {
            v$label <- seq(1, nrow(v), 1)
            for (i in 1:nrow(v)) {
                w <- v[abs(v$bin1_mid - v$bin1_mid[i]) <= GAP & abs(v$bin2_mid - v$bin2_mid[i]) <= GAP, ]
                w.min <- min(w$label)
                w.label <- sort(unique(w$label))
                for (j in 2:length(w.label)) {
                    v$label[v$label == w.label[j]] <- w.min
                }
            }
            v.rec <- sort(unique(v$label))
            v.rec <- cbind(v.rec, seq(1, length(v.rec), 1))
            for (i in 1:nrow(v)) {
                v$label[i] <- v.rec[v.rec[,1] == v$label[i], 2]
            }
            v$NegLog10P <- 0
            v$ClusterSize <- 0
            for (i in 1:nrow(v.rec)) {
                vtmp <- v[v$label == i, ]
                v$NegLog10P[v$label == i] <- sum(-log10(vtmp$p_val_reg2))
                v$ClusterSize[v$label == i] <- nrow(vtmp)
            }
            out <- rbind(out, v)
        }
        final <- rbind(final, out)
    }
    print(dim(final))
    return(final)
}

classify_peaks <- function(final) {
    if (nrow(final) == 0) return(final)
    x <- unique(final[final$label != 0, c('chr', 'label', 'NegLog10P', 'ClusterSize')])
    if (nrow(x) == 0) {
        final$ClusterType <- 'Singleton'
        return(final)
    }
    x <- x[order(x$NegLog10P), ]
    y <- sort(x$NegLog10P)
    z <- cbind(seq(1, length(y), 1), y)
    z0 <- z
    z[,1] <- z[,1] / max(z[,1])
    z[,2] <- z[,2] / max(z[,2])
    u <- z
    u[,1] <- 1/sqrt(2) * z[,1] + 1/sqrt(2) * z[,2]
    u[,2] <- -1/sqrt(2) * z[,1] + 1/sqrt(2) * z[,2]
    v <- cbind(u, seq(1, nrow(u), 1))
    RefPoint <- v[v[,2] == min(v[,2]), 3]
    RefValue <- z0[RefPoint, 2]
    final$ClusterType <- '0'
    final$ClusterType[final$label == 0] <- 'Singleton'
    final$ClusterType[final$label >= 1 & final$NegLog10P < RefValue] <- 'SharpPeak'
    final$ClusterType[final$label >= 1 & final$NegLog10P >= RefValue] <- 'BroadPeak'
    print(table(final$ClusterType))
    return(final)
}

### Main Loop
mx_combined_and = data.frame()
mx_combined_xor = data.frame()
summary_all_runs = data.frame()
singletons_names = paste(chroms, '_0', sep='')

for (r in runs) {
    name_counter = 1
    for (i in chroms) {
        print(paste('Run', r, ': regression on chromosome', i))
        mm = subset(mm_combined_and, chr == i)
        if (REG_TYPE == 'pospoisson') {
            mm = pospoisson_regression(mm, dataset_length)
        } else if (REG_TYPE == 'negbinom') {
            mm = negbinom_regression(mm, dataset_length)
        }
        write.table(mm, outf_names[name_counter], row.names=TRUE, col.names=TRUE, quote=FALSE)
        name_counter = name_counter + 1
        mx_combined_and = rbind(mx_combined_and, mm)

        mm = subset(mm_combined_xor, chr == i)
        if (REG_TYPE == 'pospoisson') {
            mm = pospoisson_regression(mm, dataset_length)
        } else if (REG_TYPE == 'negbinom') {
            mm = negbinom_regression(mm, dataset_length)
        }
        write.table(mm, outf_names[name_counter], row.names=TRUE, col.names=TRUE, quote=FALSE)
        name_counter = name_counter + 1
        mx_combined_xor = rbind(mx_combined_xor, mm)
    }

    qc_out = paste(INFDIR, SET, '.maps.qc', sep='')
    qc_val = c(sum(mm_combined_and$count, na.rm=TRUE), sum(mm_combined_xor$count, na.rm=TRUE))
    qc_label = c('AND_set', 'XOR_set')
    qc_name = c('number of sequencing pairs in AND set', 'number of sequencing pairs in XOR set')
    # Note: Original doesn’t write qc_out, so skipping unless you specify

    summary_one_run = data.frame()
    singletons = data.frame()
    for (fdr_cutoff in FDR) {
        print(paste("Finding peaks for FDR", fdr_cutoff))
        peaks_and = subset(mx_combined_and, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > fdr_cutoff)
        peaks_xor = subset(mx_combined_xor, count >= COUNT_CUTOFF & ratio2 >= RATIO_CUTOFF & -log10(fdr) > fdr_cutoff)
        if (nrow(peaks_and) > 0 && nrow(peaks_xor) > 0) {
            peaks = rbind(peaks_and, peaks_xor)
        } else if (nrow(peaks_and) > 0) {
            peaks = peaks_and
        } else if (nrow(peaks_xor) > 0) {
            peaks = peaks_xor
        } else {
            peaks = data.frame()
        }
        if (nrow(peaks) == 0) {
            print(paste('Warning: 0 bin pairs with count >= ', COUNT_CUTOFF, ' ratio >= ', RATIO_CUTOFF, ' and -log10(fdr) > ', fdr_cutoff))
        } else {
            peaks = label_peaks(peaks)
            peaks$lab = paste(peaks$chr, peaks$label, sep='_')
            peak_types = classify_peaks(peaks)
            outf_name = paste(INFDIR, SET, '.', fdr_cutoff, '.peaks', sep='')
            write.table(peak_types, outf_name, row.names=FALSE, col.names=TRUE, quote=FALSE)
            
            peak_classes = table(peaks$lab)
            n_singletons = sum(peak_classes[names(peak_classes) %in% singletons_names])
            n_singleton_chroms = length(peak_classes[names(peak_classes) %in% singletons_names])
            fraction = if (length(peak_classes) > 0) {
                n_singletons / (length(peak_classes) + n_singletons - n_singleton_chroms)
            } else {
                0
            }
            singletons_one_run = data.frame(fdr=fdr_cutoff, fraction=fraction)
            singletons = rbind(singletons, singletons_one_run)
            print(paste(fdr_cutoff, ':', fraction))
            summary_one_run = rbind(summary_one_run, do_summaries(peaks_and, peaks_xor, peaks, fraction, r))
        }
    }
    # Export all peaks (moved outside FDR loop like original intent)

    ## Original peak calling block with conditional rbind

    #if (nrow(mx_combined_and) > 0 && nrow(mx_combined_xor) > 0) {

    #    all_peaks = rbind(mx_combined_and, mx_combined_xor)

    #} else if (nrow(mx_combined_and) > 0) {

    #    all_peaks = mx_combined_and

    #} else if (nrow(mx_combined_xor) > 0) {

    #    all_peaks = mx_combined_xor

    #} else {

    #    all_peaks = data.frame()

    #}

    #if (nrow(all_peaks) > 0) {

    #    all_peaks = label_peaks(all_peaks)

    #    all_peaks$lab = paste(all_peaks$chr, all_peaks$label, sep='_')

    #    all_peak_types = classify_peaks(all_peaks)

    #    all_outf_name = paste(INFDIR, SET, '.all.peaks', sep='')

    #    fwrite(all_peak_types, all_outf_name, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

    #}
    summary_all_runs = rbind(summary_all_runs, summary_one_run)
}

summary_outf_name = paste(INFDIR, 'summary.', SET, '.txt', sep='')
write.table(summary_all_runs, summary_outf_name, row.names=FALSE, col.names=TRUE, quote=FALSE)