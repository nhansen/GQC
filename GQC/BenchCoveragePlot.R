#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/coverage_plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "year1pat")
benchname <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
genomefile <- ifelse(!( is.na(args[4])), args[4], ".")
nlocfile <- ifelse(!( is.na(args[5])), args[5], ".")
plottitle <- ifelse(!( is.na(args[6])), args[6], paste(c(benchname, " aligned coverage vs. ", genomename), sep="", collapse=""))
benchgenome <- toGRanges(genomefile)
nlocranges <- toGRanges(nlocfile)

maternalchroms <- paste0("chr", c(1:22, "X"), "_MATERNAL")
paternalchroms <- paste0("chr", c(1:22, "Y"), "_PATERNAL")
chroms <- c(maternalchroms, paternalchroms)

#benchcoveredfile <- paste(c(outputdir, "/", genomename, ".benchcovered.", benchname, ".merged.bed"), sep="", collapse="")
benchcoveredfile <- paste(c(outputdir, "/", genomename, ".benchcovered.", benchname, ".bed"), sep="", collapse="")
benchcoveredranges <- toGRanges(benchcoveredfile)
benchcovereddf <- read.table(benchcoveredfile, header=FALSE, sep="\t")

errorfile <- paste(c(outputdir, "/", genomename, ".errortype.", benchname, ".bed"), sep="", collapse="")
errorsites <- read.table(errorfile, header=FALSE, sep="\t")
names(errorsites) <- c("chrom", "start", "end", "errorname", "score", "strand", "widestart", "wideend", "rgbcolor", "errortype", "assemblyerrorname")
namefieldlengths <- sapply(seq(1, length(errorsites$chrom)), function(x) {length(strsplit(errorsites$errorname[x], split="_")[[1]])})
errorsites$ref <- sapply(seq(1,length(errorsites$chrom)), function(x) {strsplit(errorsites$errorname[x], split="_")[[1]][[namefieldlengths[x] - 1]]})
errorsites$alt <- sapply(seq(1,length(errorsites$chrom)), function(x) {strsplit(errorsites$errorname[x], split="_")[[1]][[namefieldlengths[x]]]})

# calculate error profiles for error plot:

indelerrors <- errorsites[errorsites$errortype=="CONSENSUS" & (nchar(errorsites$alt)>1 | nchar(errorsites$ref)>1 | errorsites$alt=="*" | errorsites$ref=="*"),]
indelerrorgranges <- as(indelerrors, "GRanges")
snperrors <- errorsites[errorsites$errortype=="CONSENSUS" & nchar(errorsites$alt)==1 & nchar(errorsites$ref)==1 & errorsites$alt!="*" & errorsites$ref!="*",]
snperrorgranges <- as(snperrors, "GRanges")
phasingerrors <- errorsites[errorsites$errortype=="PHASING",]
phasingerrorgranges <- as(phasingerrors, "GRanges")

plotheatmaps <- function(kp, phasingerrors, snperrors, indelerrors) {
  chromlengths <- width(ranges(benchgenome))
  names(chromlengths) <- seqnames(benchgenome)
  windowranges <- tileGenome(seqlengths=chromlengths, tilewidth=1000000)
  phasingerrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(phasingerrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=phasingerrorcounts, r0=0.1, r1=0.3, colors=c("lightblue", "blue", "black"))
  snperrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(snperrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=snperrorcounts, data.panel=1, r0=0.4, r1=0.6, colors=c("lightgreen", "darkgreen", "black"))
  indelerrorcounts <- sapply(seq(1, length(windowranges)), function(x) {sum(indelerrorgranges %over% windowranges[x])})
  kpHeatmap(kp, data=windowgranges, y=indelerrorcounts, data.panel=1, r0=0.7, r1=0.9, colors=c("white", "pink", "red"))
  
}

plotlinecurves <- function(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges) {
  kp <- kpPlotDensity(kp, data=phasingerrorgranges, data.panel=1, r0=0.1, r1=0.9, border="blue", col=NA, window.size=100000, fill=NULL, ymax=500)
  kp <- kpPlotDensity(kp, data=snperrorgranges, data.panel=1, r0=0.1, r1=0.9, col=NA, border="red", window.size=100000, ymax=500)
  kp <- kpPlotDensity(kp, data=indelerrorgranges, data.panel=1, r0=0.1, r1=0.9, col=NA, border="darkgreen", window.size=100000, ymax=500)
}

plotrectangles <- function(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges) {
  kpRect(kp, data=phasingerrorgranges, border="blue", lwd=0.1, col="blue", data.panel=1, y0=rep(0, length(phasingerrorgranges)), y1=rep(0.2, length(phasingerrorgranges)))
  # then snps:
  kpRect(kp, data=snperrorgranges, border="red", lwd=0.1, col="red", data.panel=1, y0=rep(0.2, length(snperrorgranges)), y1=rep(0.4, length(snperrorgranges)))
  # then indels:
  kpRect(kp, data=indelerrorgranges, border="darkgreen", lwd=0.1, col="darkgreen", data.panel=1, y0=rep(0.4, length(indelerrorgranges)), y1=rep(0.6, length(indelerrorgranges)))
}

plot_coverage_vs_benchmark <- function(phasingerrorgranges=NA, snperrorgranges=NA, indelerrorgranges=NA) {
  pp <- getDefaultPlotParams(plot.type=1)
  pp$ideogramheight <- 100
  pp$topmargin <- 300

  # values that were used initially:
  #pp$ideogramheight <- 400
  #pp$data1height <- 1200
  #pp$data2height <- 1000
  #pp$topmargin <- 650
  
  bordercol = 'black'
  borderlwd = 0.2
  aligncolor <- ifelse(str_detect(benchcovereddf$V1, "MAT"), "green", "blue")
  kp <- plotKaryotype(genome=benchgenome, plot.type=1, chromosomes=chroms, main=plottitle, cex=0.5, plot.params=pp)
  kpAddBaseNumbers(kp, tick.dist = 20000000, tick.len = 10, cex=0.3)
  kpRect(kp, data=benchcoveredranges, border=bordercol, lwd=borderlwd, col=aligncolor, data.panel="ideogram", y0=rep(0, length(benchcoveredranges)), y1=rep(1, length(benchcoveredranges)))
  
  kpRect(kp, data=nlocranges, border=bordercol, lwd=borderlwd, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))
  
  legend("topright", inset = c(0,  0.075), c("Paternal", "Maternal", "N's (rDNAs)", "Uncovered"), fill=c("blue", "green", "red", "grey"), title="Aligned Regions")
  
  if ((length(phasingerrorgranges)>1 || !is.na(phasingerrorgranges)) && (length(snperrorgranges)>1 || !is.na(snperrorgranges)) && (length(indelerrorgranges)>1 || !is.na(indelerrorgranges))) {
    plotlinecurves(kp, phasingerrorgranges, snperrorgranges, indelerrorgranges)
    legend("bottomright", c("Phasing", "SingleNuc", "Indel"), col=c("blue", "red", "darkgreen"), lty=c(1,1,1,1), title="Error Type")
  }
}

covgplotname <- paste(c(outputdir, "/", genomename, ".benchcovered.", benchname, ".pdf"), sep="", collapse="")
pdf(covgplotname, 8.5, 11.0)
plot_coverage_vs_benchmark()
dev.off()

errorplotname <- paste(c(outputdir, "/", genomename, ".benchcoveredwitherrors.", benchname, ".pdf"), sep="", collapse="")
pdf(errorplotname, 8.5, 11.0)
plot_coverage_vs_benchmark(phasingerrorgranges, snperrorgranges, indelerrorgranges)
dev.off()

