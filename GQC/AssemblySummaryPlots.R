
library(stringr)
args = commandArgs(trailingOnly=TRUE)

assemblyname <- ifelse(!( is.na(args[1])), args[1], "illumina2x250mat")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
idealfile <- ifelse(!( is.na(args[4])), args[4], "v1.1.non_n_seq.bed")
assemblyqv <- ifelse(!( is.na(args[5])), args[5], NA)

indelfile <- paste(c(outputdir, "/", assemblyname, ".indelerrorstats.txt"), sep="", collapse="")
subsfile <- paste(c(outputdir, "/", assemblyname, ".singlenucerrorstats.txt"), sep="", collapse="")
plotname <- paste(c(outputdir, "/", assemblyname, ".summarystats.", genomename, ".pdf"), sep="", collapse="")

#plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c("NGAx ", assemblyname, " vs ", genomename), sep="", collapse=""))
lengthtypes <- c("alignclusterlengths", "contiglengths", "scaffoldlengths")
assemblysizefiles <- sapply(lengthtypes, function(x) {file=paste(c(outputdir, "/", assemblyname, ".", x, ".txt"), sep="", collapse=""); return(file)})

mononucsitefile <- paste(c(outputdir, "/", assemblyname, ".mononucstats.txt"), sep="", collapse="")
mononucsites <- read.table(mononucsitefile, header=FALSE, sep="\t")
names(mononucsites) <- c("name", "base", "reflength", "assemblylength", "type")

consensuserrors <- mononucsites[(mononucsites$assemblylength != -1) & (mononucsites$type == "CONSENSUS"), ]
noncomplexcovered <- mononucsites[mononucsites$assemblylength != -1, ]

mononucerrorrate <- length(consensuserrors$reflength)/length(noncomplexcovered$reflength)
mononucerrorperc <- round(mononucerrorrate*100, 2)

consensuserrorcounts <- hist(consensuserrors$reflength, plot=FALSE, breaks=seq(10, 100, 1))
noncomplexcovcounts <- hist(noncomplexcovered$reflength, plot=FALSE, breaks=seq(10, 100, 1))

accrate <- 1.0 - consensuserrorcounts$counts/noncomplexcovcounts$counts

#text(20, 20, labels= paste(c("Overall error rate: ", mononucerrorperc, "%"), sep="", collapse=""))
# Make summary plot:
pdf(plotname, 11.0, 11.0)
par(mfrow=c(2,2))
assembly_ngax_plot(assemblysizefiles, assemblylabels=c("Aligned NGAx", "Contig NGx", "Scaffold NGx"), ideal=TRUE, haplotype=NA, plottitle=paste(c("Continuity stats for ", assemblyname), sep="", collapse=""))
assembly_mononucqv_plot(c(mononucsitefile), assemblylabels=c(genomename), plottitle="", plotlines=TRUE, linetype=1, errorbars=TRUE, overallerrorrate=mononucerrorperc)
assembly_compound_plot(subsfile, indelfile, readsetname, titleval=paste(c(assemblyname, " Discrepancy Counts"), sep="", collapse=""), assemblyqv=assemblyqv)
plotindellengths(indelfile, outputdir, titleval=paste(c("Indel rate by length for ", assemblyname), sep="", collapse=""))
dev.off()

