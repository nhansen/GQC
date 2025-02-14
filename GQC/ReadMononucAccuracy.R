#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/readaligns")

args = commandArgs(trailingOnly=TRUE)

readsetname <- ifelse(!( is.na(args[1])), args[1], "illumina2x250mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.0.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
plottitle <- ifelse(!( is.na(args[4])), args[4], paste(c("Homopolymer length concordance for ", readsetname, " vs ", genomename), sep="", collapse=""))

readmononuchist <- function(filename) {
  mndf <- read.table(filename, header=TRUE, sep="")
  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  mndf <- mndf[mndf$readlength != -1, ]
  mndf$reflength <- as.integer(mndf$reflength)
  mndf$readlength <- as.integer(mndf$readlength)
  mndf$numcorrect <- as.integer(mndf$numcorrect)
  mndf$numerror <- as.integer(mndf$numerror)
  
  return(mndf)
}

mndiffs <- function(filename, i) {
  mndiff <- readmononuchist(filename);
  
  if (i==0) {
    diffcount <- sum(mndiff[mndiff$readlength==mndiff$reflength, "numcorrect"])
  }
  else {
    diffcount <- sum(mndiff[mndiff$readlength-mndiff$reflength==i, "numerror"])
  }
  return(diffcount)
}

onebaseindelratio <- function(filename) {
  onebaseratio <- mndiffs(filename, 1)*1.0/mndiffs(filename, -1)
  
  return(onebaseratio)
}

diffarray <- function(filename, maxdiff=4) {
  platformdiffs <- sapply(seq(-1*maxdiff, maxdiff), function(i) {mndiffs(filename, i)})
  platformdiffs <- platformdiffs/sum(platformdiffs)
  
  return(platformdiffs)
}

makemultiplot <- function(maxdiff=4) {
  plotcolors <- c("blue", "lightblue", "lightgreen", "darkgreen", "violet", "darkviolet")
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  barplot(diffarray("sequelII_15k.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[1], cex.names=1.0, main="SequelII_15k", xlab="Difference from benchmark", ylab="Fraction of reads")
  barplot(diffarray("sequelII_20kb_20210812.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[2], cex.names=1.0, main="SequelII_20kb_20210812", xlab="Difference from benchmark")
  barplot(diffarray("hifidcv1.1.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[3], cex.names=1.0, main="HiFi/DCv1.1", xlab="Difference from benchmark")
  barplot(diffarray("hifirevio3cells.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[4], cex.names=1.0, main="HiFi Revio", xlab="Difference from benchmark")
  barplot(diffarray("illumina2x250mat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[5], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark")
  barplot(diffarray("element_avitilongmat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[6], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark")
  mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)
}

plotmononuchist <- function(file, plottitle="", maxdiff=4){
  differences <- diffarray(file, maxdiff=maxdiff)
  accrate <- as.integer(differences[maxdiff+1]*1000)/10
  barplot(differences, names=seq(-1*maxdiff,maxdiff), col="blue", cex.names=1.0, main=plottitle, xlab="Difference from benchmark", ylab="Fraction of reads", ylim=c(0,1.0))
  mtext(paste(c("Accuracy: ", as.character(accrate), "%"), sep="", collapse=""), side=3, adj=0.85, line=-4)
  indelratio <- as.integer(onebaseindelratio(file)*1000)/1000
  mtext(paste(c("1bp Ins/Dels: ", as.character(indelratio)), sep="", collapse=""), side=3, adj=0.85, line=-5 )
}

plotmononucaccuracy <- function(file, plottitle=NA, minlength=NA, maxlength=NA, color="black"){
  mncounts <- readmononuchist(file)
  if(is.na(minlength)) {
    minlength <- min(mncounts$reflength)
  }
  if(is.na(maxlength)) {
    maxlength <- max(mncounts$reflength)
  }
  if(is.na(plottitle)) {
    plottitle <- paste(c("Accuracy of mononucleotide runs in ", genomename), sep="", collapse="")
  }
  accarray <- sapply(seq(minlength, maxlength), function(x) {sum(mncounts[mncounts$reflength==x, "numcorrect"])/sum(mncounts[mncounts$reflength==x, "numerror"] + mncounts[mncounts$reflength==x, "numcorrect"])})
  #  names(mndf) <- c("reflength", "readlength", "numcorrect", "numhetallele", "numerror", "numcomplex")
  totalaccrate <- as.integer(sum(mncounts$numcorrect)*1000/sum(mncounts$numerror+mncounts$numcorrect))/10
  plot(seq(minlength, maxlength), accarray, pch=16, xlab=c("Mononucleotide run length"), ylab=c("Accuracy"), main=plottitle, ylim=c(0,1), col=color)
  if (totalaccrate < 75) {
    textheight <- 0.9
    textpos <- minlength + 3*(maxlength-minlength)/4    
  }
  else {
    #textheight <- accarray[as.integer((maxlength-minlength)/3)]/2
    textheight <- 0.2
    textpos <- minlength + (maxlength-minlength)/4
  }
  text(textpos, textheight, labels= paste(c("Overall: ", totalaccrate, "%"), sep="", collapse=""))
}

par(mfrow=c(1,1))
mononucfile <- paste(c(outputdir, "/", readsetname, ".mononuchist.txt"), sep="", collapse="")

plotname <- paste(c(outputdir, "/", readsetname, ".mononucsizehist.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotmononuchist(mononucfile, plottitle, 4)
dev.off()

plotname <- paste(c(outputdir, "/", readsetname, ".mononucaccuracy.pdf"), sep="", collapse="")
pdf(plotname, 8.5, 5.5)
plotmononucaccuracy(mononucfile, plottitle, maxlength=40)
dev.off()

comparison_plot_may_2024 <- function(maxdiff=4) {
  plotcolors <- palette.colors(n=6)
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  barplot(diffarray("element_avitilongmat.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[1], cex.names=0.7, main="Element Aviti (maternal)", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("illumina2x250mat_nonexc.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[2], cex.names=0.7, main="Illumina 2x250 (maternal)", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("hifi_revio_pbmay24.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[3], cex.names=1.0, main="HiFi Revio/Late 2023", xlab="Difference from benchmark", ylab="Fraction of reads", ylim=c(0,1))
  barplot(diffarray("ONT_Q28.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[4], cex.names=1.0, main="ONT Q28/Nov 2023", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("ONT_R10_duplex_nonexc.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[5], cex.names=1.0, main="ONT R10 Duplex", xlab="Difference from benchmark", ylim=c(0,1))
  barplot(diffarray("ONT_Q28_herrocorrected.mononuchist.txt", maxdiff=maxdiff), names=seq(-1*maxdiff, maxdiff), col=plotcolors[6], cex.names=1.0, main="Herro-corrected ONT", xlab="Difference from benchmark", ylim=c(0,1))
  mtext("Homopolymer length concordance", side = 3, line = 1, outer = TRUE)
}
accuracy_comparison_plot_may_2024 <- function() {
  plotcolors <- palette.colors(n=6)
  par(oma=c(5,5,5,5))
  par(mfrow = c(2, 3))
  plotmononucaccuracy("element_avitilongmat.mononuchist.txt", maxlength=40, color=plotcolors[1], plottitle="Element Aviti (maternal)")
  plotmononucaccuracy("illumina2x250mat_nonexc.mononuchist.txt", maxlength=40, color=plotcolors[2], plottitle="Illumina 2x250 (maternal)")
  plotmononucaccuracy("hifi_revio_pbmay24.mononuchist.txt", maxlength=40, color=plotcolors[3], plottitle="HiFi Revio/Late 2023")
  plotmononucaccuracy("ONT_Q28.mononuchist.txt", maxlength=40, color=plotcolors[4], plottitle="ONT Q28/Nov 2023")
  plotmononucaccuracy("ONT_R10_duplex_nonexc.mononuchist.txt", maxlength=40, color=plotcolors[5], plottitle="ONT R10 Duplex")
  plotmononucaccuracy("ONT_Q28_herrocorrected.mononuchist.txt", maxlength=40, color=plotcolors[6], plottitle="Herro-corrected ONT")
  mtext("Homopolymer length accuracy", side = 3, line = 1, outer = TRUE)
}

