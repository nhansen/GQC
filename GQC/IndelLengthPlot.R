#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readerrors")

library(stringr)
args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "illumina2x250mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Indel Sizes in ", readsetname, " vs ", genomename), sep="", collapse=""))

makemultiplot <- function() {
  par(mfrow=c(2,2))
  plotindellengths("sequelII_15k", outputdir, titleval="SequelII_15k")
  plotindellengths("sequelII_20kb_20210812", outputdir, titleval="SequelII_20kb_20210812")
  plotindellengths("hifidcv1.1", outputdir, titleval="HiFi/DCv1.1")
  plotindellengths("hifirevio3cell", outputdir, titleval="HiFi Revio")
  plotindellengths("illumina2x250mat", outputdir, titleval="Illumina 2x250 (maternal)")
  plotindellengths("element_avitilongmat", outputdir, titleval="Element Aviti")
}

makebog2024posterplots <- function() {
  # plots saved in "PlotsForAdamsBOGPoster2024/" called SequelDCv1.1IndelErrorRates.pdf and RevioDCv1.1.IndelErrorRates.pdf
  plotindellengths("/Users/nhansen/HG002_diploid_benchmark/plots/assemblyerrorstats/sequel_DCv1.1_hap1.indelerrorstats.txt", ".", ymax=4.0, titleval="Sequel/DCv1.1")
  plotindellengths("/Users/nhansen/HG002_diploid_benchmark/plots/assemblyerrorstats/RevioWashUDCv1.2Hap1.indelerrorstats.txt", ".", ymax=4.0, titleval="Revio/DCv1.1")
}

indelfile <- paste(c(outputdir, "/", readsetname, ".indelerrorstats.txt"), sep="", collapse="")
plotname <- paste(c(outputdir, "/", readsetname, ".indelsizestats.", genomename, ".pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotindellengths(indelfile, outputdir, titleval=plottitle)
dev.off()

