#setwd("/Users/nhansen/OneDrive/HG002_diploid_benchmark/plots/ngaplots")

#source ("./AssemblyFunctions.R")

args = commandArgs(trailingOnly=TRUE)
assemblycolors <- c("#44AA99", "#332288", "#882255", "#888888", "#DDCC77", "#661100", "#6699CC")
methodcolors <- c("#DDCC77", "#661100", "#6699CC")

assemblyname <- ifelse(!( is.na(args[1])), args[1], "lc24_ul_herro_ulk")
genomename <- ifelse(!( is.na(args[2])), args[2], "v1.1")
outputdir <- ifelse(!( is.na(args[3])), args[3], ".")
idealfile <- ifelse(!( is.na(args[4])), args[4], "v1.1.non_n_seq.bed")
plottitle <- ifelse(!( is.na(args[5])), args[5], paste(c("NGAx ", assemblyname, " vs ", genomename), sep="", collapse=""))
lengthtypes <- c("alignclusterlengths", "contiglengths", "scaffoldlengths")
assemblysizefiles <- sapply(lengthtypes, function(x) {file=paste(c(outputdir, "/", assemblyname, ".", x, ".txt"), sep="", collapse=""); return(file)})

# Make NGAx plot:
plotname <- paste(c(outputdir, "/", assemblyname, ".continuitystats.", genomename, ".pdf"), sep="", collapse="")
pdf(plotname)
assembly_ngax_plot(assemblysizefiles, assemblylabels=c("Aligned NGAx", "Contig NGx", "Scaffold NGx"), ideal=TRUE, haplotype=NA, plottitle=paste(c("Continuity stats for ", assemblyname), sep="", col=""))
dev.off()

