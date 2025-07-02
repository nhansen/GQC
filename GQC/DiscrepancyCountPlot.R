#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readerrors")

library(stringr)
args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "illumina2x250mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Discrepancy counts in ", readsetname, " vs ", genomename), sep="", collapse=""))

indelfile <- paste(c(outputdir, "/", readsetname, ".indelerrorstats.txt"), sep="", collapse="")
subsfile <- paste(c(outputdir, "/", readsetname, ".singlenucerrorstats.txt"), sep="", collapse="")
plotname <- paste(c(outputdir, "/", readsetname, ".discrepancystats.", genomename, ".pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
assembly_compound_plot(subsfile, indelfile, readsetname, titleval=plottitle)
dev.off()

