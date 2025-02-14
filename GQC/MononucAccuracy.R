#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/mononucplots")

library("Hmisc")

args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "ash1")
benchname <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
resourcedir <- ifelse(! ( is.na(args[4])), args[4], ".")
plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c(genomename, " accuracy of mononucleotide runs"), sep="", collapse=""))

mononucsitefile <- paste(c(outputdir, "/", genomename, ".mononucstats.txt"), sep="", collapse="")
mononucsites <- read.table(mononucsitefile, header=FALSE, sep="\t")
names(mononucsites) <- c("name", "base", "reflength", "assemblylength", "type")

consensuserrors <- mononucsites[(mononucsites$assemblylength != -1) & (mononucsites$type == "CONSENSUS"), ]
noncomplexcovered <- mononucsites[mononucsites$assemblylength != -1, ]

mononucerrorrate <- length(consensuserrors$reflength)/length(noncomplexcovered$reflength)
mononucerrorperc <- round(mononucerrorrate*100, 2)

consensuserrorcounts <- hist(consensuserrors$reflength, plot=FALSE, breaks=seq(10, 100, 1))
noncomplexcovcounts <- hist(noncomplexcovered$reflength, plot=FALSE, breaks=seq(10, 100, 1))

accrate <- 1.0 - consensuserrorcounts$counts/noncomplexcovcounts$counts

plotname <- paste(c(outputdir, "/", genomename, ".mononuc_accuracy.", benchname, ".pdf"), sep="", collapse="")
pdf(plotname, 7.08661, 6.69291)

plot(consensuserrorcounts$mids, accrate, xlim=c(10,40), pch=16, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=paste(c("Accuracy of mononucleotide runs in ", genomename), sep="", collapse=""))
text(20, 0.5, labels= paste(c("Overall error rate: ", mononucerrorperc, "%"), sep="", collapse=""))

dev.off()

