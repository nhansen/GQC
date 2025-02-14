#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/qualscores")

args = commandArgs(trailingOnly=TRUE)

genomename <- ifelse(!( is.na(args[1])), args[1], "lc24_medaka_hap1")
benchname <- ifelse(!( is.na(args[1])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
resourcedir <- ifelse(! ( is.na(args[4])), args[4], ".")
title <- ifelse(!( is.na(args[5])), args[5], paste(c(genomename, " accuracy of quality scores"), sep="", collapse=""))

qverrorsfile <- paste(c(outputdir, "/", genomename, ".qvstats.txt"), sep="", collapse="")

qvalues <- function(errorlist, totallist) {
  qscores <- sapply(seq(1, length(errorlist)), function(x) {return(ifelse(totallist[x]==0, NA, as.integer(-10.0*log10((errorlist[x]+1)/(totallist[x]+1)))))})
  return(qscores)
}

get_qv_counts <- function(file) {
  qvcounts <- read.table(file, header=FALSE, sep="\t")
  names(qvcounts) <- c("QVReported", "SNVErrors", "IndelErrors", "TotalBases")

  qvcounts$SNVObsQV <- qvalues(qvcounts$SNVErrors, qvcounts$TotalBases)
  qvcounts$IndelObsQV <- qvalues(qvcounts$IndelErrors, qvcounts$TotalBases)
  qvcounts$TotalObsQV <- qvalues((qvcounts$SNVErrors+qvcounts$IndelErrors), qvcounts$TotalBases)
  
  return(qvcounts)
}

qvcounts <- get_qv_counts(qverrorsfile)
condenseqvcounts <- function(qvcounts, windowsize=5) {
  numwindows <- as.integer(max(qvcounts$QVReported)/windowsize) + 1
  condensedcounts$SNVErrors <- sapply()
  return(numwindows)
}

makeqvplot <- function(qvcounts, cexfactor=0.05, plottitle="") {
  maxreported <- max(qvcounts[qvcounts$TotalBases>0, "QVReported"])
  plot(c(0,maxreported), c(0,maxreported), type="l", lty=2, main=plottitle, xlab="Reported QV", ylab="Observed QV")
  points(qvcounts$QVReported, qvcounts$TotalObsQV, pch=16, cex=cexfactor*log(qvcounts$TotalBases), col="black")
  points(qvcounts$QVReported, qvcounts$SNVObsQV, pch=16, cex=cexfactor*log(qvcounts$TotalBases), col="darkgreen")
  points(qvcounts$QVReported, qvcounts$IndelObsQV, pch=16, cex=cexfactor*log(qvcounts$TotalBases), col="blue")  #with(qvcounts, symbols(x=QVReported, y=TotalObsQV, circles=log(TotalBases), add=TRUE, inches=0.05, ann=F, bg="black", fg=NULL))
  legend(5, maxreported-5, c("Total QV", "SNV QV", "Indel QV"), pch=16, col=c("black", "darkgreen", "blue"))
  sizelegendsizes <- (c(2000, 200000, 20000000, 200000000, 2000000000))
  sizecexvals <- sapply(sizelegendsizes, function(x) {return(cexfactor*log(x))})
  legend(50, 20, as.character(sizelegendsizes), pch=16, col="black", pt.cex=sizecexvals, title="Total Bases")
}

plotname <- paste(c(outputdir, "/", genomename, ".qualvalue_accuracy_vs_", benchname, ".pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
makeqvplot(qvcounts, plottitle=title)
dev.off()

