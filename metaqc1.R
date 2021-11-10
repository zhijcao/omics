

#library(filesstrings)

## Use socket based parallel processing on Windows systems
## if (.Platform$OS.type == "unix") {
##     register(bpstart(MulticoreParam(2)))
## } else {
##     register(bpstart(SnowParam(2)))
## }
register(SerialParam()) 

## ----load-libs-pheno, message = FALSE--------------------------------------
library(xcms)
#library(faahKO)
library(RColorBrewer)
library(pander)
library(magrittr)
library(dplyr)
library(tidyr)

library(devtools)
library(installr)
devtools::install_github("sneumann/xcms")
## Get the full path to the CDF files
rm(list=ls())

help("installr")
installr()
library(filesstrings)

serum <- read.csv("C:/zhijuncao/metabolomics/serum samples.csv", stringsAsFactors = FALSE)
urine <- read.csv("C:/zhijuncao/metabolomics/urine samples.csv", stringsAsFactors = FALSE)
nist <- read.csv("C:/zhijuncao/metabolomics/nist samples.csv", stringsAsFactors = FALSE)

fldpath <- "C:/zhijuncao/metabolomics/QC Converted/"
fldpath1 <- "C:/zhijuncao/metabolomics/NIST Converted data/"

fldpath2 <- "C:/zhijuncao/metabolomics/qc/"

serum[,"cdfs"] <-paste(fldpath, serum$Name,"01.CDF", sep="")
urine[,"cdfs"]  <-paste(fldpath, urine$Name,"01.CDF", sep="")
nist[,"cdfs"] <- paste(fldpath1, nist$Name,"01.CDF", sep="")

serum[,"Name.raw"] <-paste(fldpath2, serum$Name,".raw", sep="")
urine[,"Name.raw"]  <-paste(fldpath2, urine$Name,".raw", sep="")
nist[,"Name.raw"] <- paste(fldpath2, nist$Name,".raw", sep="")

serump <- serum %>% filter(PosNeg=="serum")
serumn <- serum %>% filter(PosNeg=="serum_neg")

names(serum)
# for (i in (2:31)){
# file.copy(serump$Name.raw[i], "C:/zhijuncao/metabolomics/qc/serumpos/",
#           recursive = TRUE)
# }
# 
# for (i in (1:31)){
#   file.copy(serumn$Name.raw[i], "C:/zhijuncao/metabolomics/qc/serumneg/",
#             recursive = TRUE)
# }


urinep <- urine %>% filter(PosNeg=="Full scan pos2fun")
urinen <- urine %>% filter(PosNeg=="Full scan neg2fun") %>% filter(Name!="NEG_Q-tof-10-0226")
# urinen$Name.raw
# for (i in (1:34)){
#   file.copy(urinep$Name.raw[i], "C:/zhijuncao/metabolomics/qc/urinepos/",
#             recursive = TRUE)
# }
# 
# for (i in (1:27)){
#   file.copy(urinen$Name.raw[i], "C:/zhijuncao/metabolomics/qc/urineneg/",
#             recursive = TRUE)
# }


nistp <- nist %>% filter(PosNeg=="serum")
nistp1 <- nist %>% filter(PosNeg=="serum") %>% filter(Project!="CCl4_P") #remove CCl4_P, it is serum not plasma
nistn <- nist %>% filter(PosNeg=="serum_N")
# nistp$Name.raw
# nistn$Name.raw
# for (i in (1:17)){
#   file.copy(nistp$Name.raw[i], "C:/zhijuncao/metabolomics/qc/nistpos/",
#             recursive = TRUE)
# }
# 
# for (i in (1:18)){
#   file.copy(nistn$Name.raw[i], "C:/zhijuncao/metabolomics/qc/nistneg/",
#             recursive = TRUE)
# }







## Create a phenodata data.frame

serumppd <- data.frame(sample_name=basename(serump$cdfs),
                   sample_group=serump$Project, stringsAsFactors=FALSE)

serumnpd <- data.frame(sample_name=basename(serumn$cdfs),
                       sample_group=serumn$Project, stringsAsFactors=FALSE)

urineppd <- data.frame(sample_name=basename(urinep$cdfs),
                       sample_group=urinep$Project, stringsAsFactors=FALSE)

urinenpd <- data.frame(sample_name=basename(urinen$cdfs),
                       sample_group=urinen$Project, stringsAsFactors=FALSE)

nistppd <- data.frame(sample_name=basename(nistp$cdfs),
                       sample_group=nistp$Project, stringsAsFactors=FALSE)

nistp1pd <- data.frame(sample_name=basename(nistp1$cdfs),
                      sample_group=nistp1$Project, stringsAsFactors=FALSE)

nistnpd <- data.frame(sample_name=basename(nistn$cdfs),
                     sample_group=nistn$Project, stringsAsFactors=FALSE)

## ----load-with-msnbase, message = FALSE------------------------------------

serump_raw <- readMSData(files =serump$cdfs, pdata = new("NAnnotatedDataFrame", serumppd),
                         mode = "onDisk") 
saveRDS(serump_raw, "C:/zhijuncao/metabolomics/lcmsre/serump_raw.rds")

serumn_raw <- readMSData(files =serumn$cdfs, pdata = new("NAnnotatedDataFrame", serumnpd),
                         mode = "onDisk") 
saveRDS(serumn_raw, "C:/zhijuncao/metabolomics/lcmsre/serumn_raw.rds")

urinep_raw <- readMSData(files =urinep$cdfs, pdata = new("NAnnotatedDataFrame", urineppd),
                         mode = "onDisk") 
saveRDS(urinep_raw, "C:/zhijuncao/metabolomics/lcmsre/urinep_raw.rds")

urinen_raw <- readMSData(files =urinen$cdfs, pdata = new("NAnnotatedDataFrame", urinenpd),
                         mode = "onDisk") 
saveRDS(urinen_raw, "C:/zhijuncao/metabolomics/lcmsre1/urinen_raw.rds")
saveRDS(urinen_raw, "C:/zhijuncao/metabolomics/lcmsre1/urinen_raw_rm24.rds")

nistp_raw <- readMSData(files = nistp$cdfs, pdata = new("NAnnotatedDataFrame", nistppd),
                       mode = "onDisk") 
saveRDS(nistp_raw, "C:/zhijuncao/metabolomics/lcmsre1/nistp_raw.rds")
nistp_raw <- readRDS("C:/zhijuncao/metabolomics/lcmsre1/nistp_raw.rds")

nistp1_raw <- readMSData(files = nistp1$cdfs, pdata = new("NAnnotatedDataFrame", nistp1pd),
                        mode = "onDisk") 
saveRDS(nistp1_raw, "C:/zhijuncao/metabolomics/lcmsre1/nistp1_raw.rds")
nistp1_raw <- readRDS("C:/zhijuncao/metabolomics/lcmsre1/nistp1_raw.rds")

nistn_raw <- readMSData(files = nistn$cdfs, pdata = new("NAnnotatedDataFrame", nistnpd),
                       mode = "onDisk") 
saveRDS(nistn_raw, "C:/zhijuncao/metabolomics/lcmsre1/nistn_raw.rds")

group_colors <- c("red", "black","blue","cyan","gray", "pink", "yellow", "green")


##########################################################################
#############     start analysis
##########################################################################


######serum_pos
pd <- serumppd
raw <- serump_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/serum_pos"
temname <- "serum_pos"

######serum_neg
pd <- serumnpd
raw <- serumn_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/serum_neg"
temname <- "serum_neg"


######urine_pos
pd <- urineppd
raw <- urinep_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/urine_pos"
temname <- "urine_pos"

######urine_neg
pd <- urinenpd
raw <- urinen_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/urine_neg_rm24"
temname <- "urine_neg_rm24"

######nist_pos
pd <- nistppd
raw <- nistp_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/nist_pos"
temname <- "nist_pos"

### remove project CCL4_P (serum sample, not plasma)
pd <- nistp1pd
raw <- nistp1_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre1/nist1_pos"
temname <- "nist1_pos"


#######nist_neg
pd <- nistnpd
raw <- nistn_raw

pdfname <- "C:/zhijuncao/metabolomics/lcmsre/nist_neg"
temname <- "nist_neg"



#############

names(group_colors)<- unique(pd$sample_group)

## ----data-inspection-bpc, message = FALSE, fig.align = "center", fig.width = 12, fig.height = 6----
## Get the base peak chromatograms. This reads data from the files.
raw_bp <- chromatogram(raw, aggregationFun = "max")
str(raw_bp)
## Plot all chromatograms.
#str(nist_bp)

pdf(paste(pdfname,"_raw_basepeak.pdf", sep=""), width = 10, height = 6)
plot(raw_bp, col = group_colors[raw$sample_group])
legend("topright", xpd = TRUE,legend=unique(pd$sample_group), col=group_colors[unique(pd$sample_group)],  lty=1)
dev.off()


## ----data-inspection-tic-boxplot, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 4, fig.cap = "Distribution of total ion currents per file."----
## Get the total ion current by file
raw_tc <- split(tic(raw), f = fromFile(raw))

pdf(paste(pdfname,"_raw_tc.pdf", sep=""), width = 10, height = 6)
boxplot(raw_tc, col = group_colors[raw$sample_group],
        ylab = "intensity", main = "Total ion current")
legend("topright", xpd = TRUE,legend=unique(pd$sample_group), col=group_colors[unique(pd$sample_group)], pch=16)
dev.off()



## ----peak-detection-centwave, message = FALSE, results = "hide"------------
cwp <- CentWaveParam(ppm = 25, peakwidth = c(5, 20), snthresh = 10,
                     prefilter = c(0, 0), mzCenterFun = "wMean", integrate = 1L,
                     mzdiff = -0.001, fitgauss = FALSE, noise = 10, verboseColumns = FALSE,
                     roiList = list(), firstBaselineCheck = TRUE, roiScales = numeric())

xraw<- findChromPeaks(raw, param = cwp) 


## ----peak-detection-peaks-per-sample, message = FALSE, results = "asis"----
summary_fun <- function(z) {
  c(peak_count = nrow(z), rt = quantile(z[, "rtmax"] - z[, "rtmin"]))
}
peakinformation <- lapply(split.data.frame(chromPeaks(xraw),
                             f = chromPeaks(xraw)[, "sample"]),
            FUN = summary_fun)


peakinformation <- do.call(rbind, peakinformation)
rownames(peakinformation) <- basename(fileNames(xraw))
pandoc.table(peakinformation,
             caption = paste0("Summary statistics on identified chromatographic",
                              " peaks. Shown are number of identified peaks per",
                              " sample and widths/duration of chromatographic ",
                              "peaks.")) 
write.csv(peakinformation,paste(pdfname,"_peakinformation.csv", sep=""))


## ----peak-detection-chrom-peak-image, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Frequency of identified chromatographic peaks along the retention time axis. The frequency is color coded with higher frequency being represented by yellow-white. Each line shows the peak frequency for one file."----
pdf(paste(pdfname,"_xraw_image.pdf", sep=""), width = 6, height = 6)
plotChromPeakImage(xraw)
dev.off()



## ----peak-detection-chrom-peak-intensity-boxplot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 8, fig.cap = "Peak intensity distribution per sample."----
## Extract a list of per-sample peak intensities (in log2 scale)
ints <- split(log2(chromPeaks(xraw)[, "into"]),
              f = chromPeaks(xraw)[, "sample"])

pdf(paste(pdfname,"_xraw_peaks_boxplot.pdf", sep=""), width = 10, height = 6)
boxplot(ints, varwidth = TRUE, col = group_colors[xraw$sample_group],
        ylab = expression(log[2]~intensity), main = "Peak intensities")
grid(nx = NA, ny = NULL) 
dev.off()

#############################################################################
## ----alignment-obiwarp, message = FALSE, results = "hide"------------------
#############################################################################
xraw <- adjustRtime(xraw, param = ObiwarpParam(binSize = 0.6)) 


## ----alignment-obiwarp-plot, message = FALSE, fig.align = "center", fig.width = 10, fig.height = 10, fig.cap = "Obiwarp aligned data. Base peak chromatogram after alignment (top) and difference between adjusted and raw retention times along the retention time axis (bottom)."----
## Get the base peak chromatograms.
xraw_bp_adj <- chromatogram(xraw, aggregationFun = "max")
str(xraw_bp_adj)

pdf(paste(pdfname,"_raw_bp_xraw_adj_bp.pdf", sep=""), width = 10, height = 9)
par(mfrow = c(3, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(raw_bp, col = group_colors[raw$sample_group])
legend("topright", xpd = TRUE,legend=unique(pd$sample_group), col=group_colors[unique(pd$sample_group)],  lty=1)

plot(xraw_bp_adj , col = group_colors[xraw_bp_adj$sample_group], xlab="adj retention time")
## Plot also the difference of adjusted to raw retention time.
plotAdjustedRtime(xraw , col = group_colors[xraw$sample_group]) 
dev.off()
#phenylalanine:120.08, 166.08
mzr <- c(119.8,120.18)
mzr <- c(165.98,166.18)
rtr <- c(120,150)
chr_raw <- chromatogram(raw, mz = mzr, rt = rtr, aggregationFun = "sum")
chr_xraw <- chromatogram(xraw, mz = mzr, rt = rtr, aggregationFun = "sum")

pdf(paste(pdfname,"_mz_RT.pdf", sep=""), width = 4, height = 3.5)
plotChromPeaks(xraw, file = 1)
dev.off()

pdf(paste(pdfname,"_raw_bp_xraw_adj_bp_phenylalanine166.pdf", sep=""), width = 10, height = 9)
par(mfrow = c(3, 1), mar = c(4.5, 4.2, 1, 0.5))
plot(chr_raw, col = group_colors[raw$sample_group])
legend("topright", xpd = TRUE,legend=unique(pd$sample_group), col=group_colors[unique(pd$sample_group)],  lty=1)

plot(chr_xraw , col = group_colors[xraw_bp_adj$sample_group], xlab="adj retention time")
## Plot also the difference of adjusted to raw retention time.
plotAdjustedRtime(xraw , col = group_colors[xraw$sample_group]) 
dev.off()




## ----alignment-peak-groups, message = FALSE--------------------------------
## Correspondence: group peaks across samples.
#################################################

pdp <- PeakDensityParam(sampleGroups = xraw$sample_group,
                        minFraction = 0.8, binSize = 0.015, bw = 2)
xraw <- groupChromPeaks(xraw, param = pdp)



## ----fill-chrom-peaks, message = FALSE-------------------------------------
## Filling missing peaks using default settings. Alternatively we could
## pass a FillChromPeaksParam object to the method.
xraw <- fillChromPeaks(xraw)


## ----fill-chrom-peaks-compare, message = FALSE-----------------------------
## Missing values before filling in peaks
apply(featureValues(xraw, filled = FALSE), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))

## Missing values after filling in peaks
apply(featureValues(xraw), MARGIN = 2,
      FUN = function(z) sum(is.na(z)))

## ----final-pca, message = FALSE, fig.align = "center", fig.width = 8, fig.height = 8, fig.cap = "PCA for the faahKO data set, un-normalized intensities."----
## Extract the features and log2 transform them
ft_ints <- featureValues(xraw, value = "into")
ftd <- featureDefinitions(xraw)
ft <- as.data.frame(cbind(ftd, ft_ints))
ft$peakidx <- as.character(ft$peakidx)
write.csv(ft, paste(pdfname,"_feature_table.csv"))

ft$peakidx <- as.character(ft$peakidx)
saveRDS(cbind(ftd, ft_ints), paste(pdfname,"_feature_table.rds"))


## Perform the PCA omitting all features with an NA in any of the
## samples. Also, the intensities are mean centered.
pc <- prcomp(t(log2(na.omit(ft_ints))), center = TRUE)

## Plot the PCA
cols <- group_colors[xraw$sample_group]
pcSummary <- summary(pc)

pdf(paste(pdfname,"_xraw_adj_goup_PCA.pdf", sep=""), width = 6, height = 6)
plot(pc$x[, 1], pc$x[,2], pch = 21, main = temname, 
     xlab = paste0("PC1: ", format(pcSummary$importance[2, 1] * 100,
                                   digits = 3), " % variance"),
     ylab = paste0("PC2: ", format(pcSummary$importance[2, 2] * 100,
                                   digits = 3), " % variance"),
     col = "darkgrey", bg = cols, cex = 2)
grid()
legend("topright", xpd = TRUE,legend=unique(pd$sample_group), col=group_colors[unique(pd$sample_group)], pch=16)
# text(pc$x[, 1], pc$x[,2], labels = xraw$sample_name, col = "darkgrey",
#      pos = 3, cex = 2)
dev.off()


####save processed chromatogram file

saveRDS(xraw, paste(pdfname,"_xraw_rm24.rds", sep=""))

#####
nist_pos_feature_table <- as.data.frame(readRDS("nist_pos_feature_table.rds"))
nist_pos_feature_table$peakidx <- as.character(nist_pos_feature_table$peakidx)
write.csv(nist_pos_feature_table, "nist_pos_feature_table.csv")

nist_pos_feature_table <- read.csv("C:/zhijuncao/metabolomics/lcmsre/nist_pos_feature_table.csv", check.names = FALSE)
nist_pos_ft <- nist_pos_feature_table %>% gather(key=Replicate, value=Total.Area, -c(metaID:peakidx))
nist_pos_ft1 <- left_join(nist_pos_ft, nist, by=c("Replicate"="Name"))
write.csv(nist_pos_ft1, "C:/zhijuncao/metabolomics/lcmsre/nist_pos_feature_table_long.csv")


nist_neg_feature_table <- as.data.frame(readRDS("nist_neg_feature_table.rds"))
nist_neg_feature_table$peakidx <- as.character(nist_neg_feature_table$peakidx)
write.csv(nist_neg_feature_table, "nist_neg_feature_table.csv")


serum_pos_feature_table <- as.data.frame(readRDS("serum_pos_feature_table.rds"))
serum_pos_feature_table$peakidx <- as.character(serum_pos_feature_table$peakidx)
write.csv(serum_pos_feature_table, "serum_pos_feature_table.csv")

serum_pos_feature_table <- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_pos_feature_table.csv", check.names = FALSE)
serum_pos_ft <- serum_pos_feature_table %>% gather(key=Replicate, value=Total.Area, -c(metaID:peakidx))
serum_pos_ft1 <- left_join(serum_pos_ft, serum, by=c("Replicate"="Name"))
write.csv(serum_pos_ft1, "C:/zhijuncao/metabolomics/lcmsre/serum_pos_feature_table_long.csv")


serum_neg_feature_table <- readRDS("serum_neg_feature_table.rds")
serum_neg_feature_table <- as.data.frame(readRDS("serum_neg_feature_table.rds"))

serum_neg_feature_table$peakidx <- as.character(serum_neg_feature_table$peakidx)
write.csv(serum_neg_feature_table, "serum_neg_feature_table.csv")

names(serum_neg_feature_table)
serum_neg_feature_table <- read.csv("C:/zhijuncao/metabolomics/lcmsre/serum_neg_feature_table.csv", check.names = FALSE)
serum_neg_ft <- serum_neg_feature_table %>% gather(key=Replicate, value=Total.Area, -c(metaID:peakidx))
serum_neg_ft1 <- left_join(serum_neg_ft, serum, by=c("Replicate"="Name"))
write.csv(serum_neg_ft1, "C:/zhijuncao/metabolomics/lcmsre/serum_neg_feature_table_long.csv")

urine_pos_feature_table <- as.data.frame(readRDS("urine_pos_feature_table.rds"))
urine_pos_feature_table$peakidx <- as.character(urine_pos_feature_table$peakidx)
write.csv(urine_pos_feature_table, "urine_pos_feature_table.csv")

urine_neg_feature_table <- as.data.frame(readRDS("urine_neg_feature_table.rds"))
urine_neg_feature_table$peakidx <- as.character(urine_neg_feature_table$peakidx)
write.csv(urine_neg_feature_table, "urine_neg_feature_table.csv")



ft <- nist_pos_feature_table
ft1 <- as.data.frame(t(ft))

head(ft1, 2)
str(ft1)



######





















## ----processhistory, message = FALSE---------------------------------------
processHistory(xdata) 
processHistory(xnistp)

## ----processhistory-select, message = FALSE--------------------------------
ph <- processHistory(xdata, type = "Retention time correction")

ph 

## ----processhistory-param, message = FALSE---------------------------------
## Access the parameter
processParam(ph[[1]])


## ----subset-filterFile, message = FALSE------------------------------------
subs <- filterFile(xdata, file = c(2, 4))

## Do we have identified chromatographic peaks?
hasChromPeaks(subs) 

## ----subset-filterFile-2, message = FALSE----------------------------------
## Do we still have features?
hasFeatures(subs)


## Do we still have adjusted retention times?
hasAdjustedRtime(subs) 

## ----subset-filterFile-3, message = FALSE----------------------------------
subs <- filterFile(xdata, keepAdjustedRtime = TRUE)

hasAdjustedRtime(subs) 

## ----subset-filterRt, message = FALSE--------------------------------------
subs <- filterRt(xdata, rt = c(3000, 3500))


range(rtime(subs)) 

## ----subset-filterRt-2, message = FALSE------------------------------------
hasAdjustedRtime(subs) 

## ----subset-filterRt-3, message = FALSE------------------------------------
hasChromPeaks(subs)

range(chromPeaks(subs)[, "rt"]) 

## ----subset-bracket, message = FALSE, warning = FALSE----------------------
## Extract all data from the 3rd file.
one_file <- filterFile(xdata, file = 3)

one_file_2 <- xdata[fromFile(xdata) == 3]

## Is the content the same?
all.equal(spectra(one_file), spectra(one_file_2)) 

## ----subset-bracket-peaks, message = FALSE---------------------------------
## Subsetting with filterFile preserves chromatographic peaks
head(chromPeaks(one_file))

## Subsetting with [ not
head(chromPeaks(one_file_2)) 

## ----subset-bracket-keepAdjustedRtime, message = FALSE, warnings = FALSE----
subs <- xdata[20:30, keepAdjustedRtime = TRUE]

hasAdjustedRtime(subs)

## Access adjusted retention times:
rtime(subs)

## Access raw retention times:
rtime(subs, adjusted = FALSE) 

## ----subset-double-bracket, message = FALSE--------------------------------
## Extract a single spectrum
xdata[[14]] 

## ----subset-split, message = FALSE-----------------------------------------
x_list <- split(xdata, f = fromFile(xdata), keepAdjustedRtime = TRUE)

lengths(x_list)

lapply(x_list, hasAdjustedRtime) 

## ----subset-split-by-file, message = FALSE---------------------------------
xdata_by_file <- splitByFile(xdata, f = factor(1:length(fileNames(xdata))))

lapply(xdata_by_file, hasChromPeaks) 

## ----multicore, message = FALSE, eval = FALSE------------------------------
#  register(bpstart(MulticoreParam(2)))

## ----snow, message = FALSE, eval = FALSE-----------------------------------

#  register(bpstart(SnowParam(2)))

