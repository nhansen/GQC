#setwd("/Users/nhansen/HG002_diploid_benchmark/plots")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#if (!requireNamespace("karyoploteR", quietly = TRUE))
  #BiocManager::install("karyoploteR")

library(stringr)
library(karyoploteR)
args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "year1pat")
outputdir <- ifelse(!( is.na(args[2])), args[2], ".")
resourcedir <- ifelse(! ( is.na(args[3])), args[3], ".")
benchname <- ifelse(!( is.na(args[4])), args[4], "v1.1")
plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c(genomename, " aligned coverage vs. v1.1"), sep="", collapse=""))
genomefile <- paste(c(outputdir, "/genome.", genomename, ".bed"), sep="", collapse="")
scaffolds <- toGRanges(genomefile)
scaffolds <- sort(scaffolds)
scaffolddf <- read.table(genomefile, sep="\t", comment.char="", header=FALSE)
seqlengths(scaffolds) <- scaffolddf[,3]
bigscaffolds <- sort(seqlengths(scaffolds), decreasing=TRUE)
chroms <- names(bigscaffolds)
if (length(chroms) > 40){
  chroms <- chroms[1:40]
}

#matcoveredfile <- paste(c(outputdir, "/testmatcovered.", genomename, ".merged.bed"), sep="", collapse="")
#patcoveredfile <- paste(c(outputdir, "/testpatcovered.", genomename, ".merged.bed"), sep="", collapse="")
matcoveredfile <- paste(c(outputdir, "/testmatcovered.", genomename, ".bed"), sep="", collapse="")
patcoveredfile <- paste(c(outputdir, "/testpatcovered.", genomename, ".bed"), sep="", collapse="")
matcoveredranges <- toGRanges(matcoveredfile)
patcoveredranges <- toGRanges(patcoveredfile)

nlocfile <- paste(c(outputdir, "/nlocs.", genomename, ".bed"), sep="", collapse="")
nlocranges <- toGRanges(nlocfile)

pp <- getDefaultPlotParams(plot.type=1)
pp$ideogramheight <- 100
pp$topmargin <- 300

plotname <- paste(c(outputdir, "/", genomename, ".testcovered.", benchname, ".pdf"), sep="", collapse="")
pdf(plotname, 8.5, 11.0)
#kp <- plotKaryotype(genome=scaffolds,chromosomes=paste0("chr", c(1:22, "X", "Y")), plot.type=1, main=plottitle)
kp <- plotKaryotype(genome=scaffolds, plot.type=1, chromosomes=chroms, main=plottitle, cex=0.5, plot.params=pp)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red")
if (length(nlocranges) > 0)
    kpRect(kp, data=nlocranges, col="red", data.panel="ideogram", y0=rep(0, length(nlocranges)), y1=rep(1, length(nlocranges)))
if (length(patcoveredranges) > 0)
    kpRect(kp, data=patcoveredranges, col="blue", data.panel="ideogram", y0=rep(0, length(patcoveredranges)), y1=rep(1, length(patcoveredranges)))
if (length(matcoveredranges) > 0)
    kpRect(kp, data=matcoveredranges, col="green", data.panel="ideogram", y0=rep(0, length(matcoveredranges)), y1=rep(1, length(matcoveredranges)))

legend("bottomright", c("Paternal", "Maternal", "Ns"), fill=c("blue", "green", "red"))

dev.off()
