library(stringr)

#setwd("/Users/nhansen/HG002_diploid_benchmark/plots/align_plots")
args = commandArgs(trailingOnly=TRUE)

#pardefault <- par()

chromfile <- ifelse(!( is.na(args[1])), args[1], "clustered_aligns.chr1_MATERNAL.clusters.bed")
genomename <- ifelse(!( is.na(args[2])), args[2], "year1pat")
benchname <- ifelse(!( is.na(args[3])), args[3], "v1.0.1")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")
resourcedir <- ifelse(! ( is.na(args[5])), args[5], ".")
chromlength <- ifelse(! ( is.na(args[6])), as.integer(args[6]), NA)

readaligns <- function(chromfile) {
  aligns <- read.table(chromfile, sep="\t", header=FALSE)
  names(aligns) <- c("chrom", "start", "end", "alignname")

  aligns$start <- as.numeric(aligns$start)/1000000
  aligns$end <- as.numeric(aligns$end)/1000000
  aligns$query <- sapply(seq(1, length(aligns$alignname)), function(i) {strsplit(aligns$alignname, split="_")[[i]][1]})
  aligns$querystart <- sapply(seq(1, length(aligns$alignname)), function(i) {as.numeric(strsplit(aligns$alignname, split="_")[[i]][2])/1000000})
  aligns$queryend <- sapply(seq(1, length(aligns$alignname)), function(i) {as.numeric(strsplit(aligns$alignname, split="_")[[i]][3])/1000000})
  aligns$cluster <- sapply(seq(1, length(aligns$alignname)), function(i) {strsplit(aligns$alignname, split="_")[[i]][4]})

  chromorderedaligns <- aligns[order(aligns$start, aligns$end), ]
  
  return(chromorderedaligns)
}

minpos <- function(aligns, contigname) {
  startsends <- c(aligns[aligns$query==contigname, "querystart"], aligns[aligns$query==contigname, "queryend"])
  return(min(startsends))
}

maxpos <- function(aligns, contigname) {
  startsends <- c(aligns[aligns$query==contigname, "querystart"], aligns[aligns$query==contigname, "queryend"])
  return(max(startsends))
}

contiginfo <- function(aligns) {
  
  ctginfo <- data.frame("contigname"=unique(aligns$query))
  ctginfo$minpos <- sapply(ctginfo$contigname, function(x) {minpos(aligns, x)})
  ctginfo$maxpos <- sapply(ctginfo$contigname, function(x) {maxpos(aligns, x)})
  ctginfo$basescovered <- ctginfo$maxpos - ctginfo$minpos
  
  return(ctginfo)
}


orderedcontigs <- function(aligns) {
  ctginfo <- contiginfo(aligns) 
  orderedctginfo <- ctginfo[order(ctginfo$basescovered, decreasing = TRUE), ]
  
  return(orderedctginfo)
}

plotaligns <- function(aligns, index=1, suppressaxis=FALSE, suppresstitle=FALSE, suppressxunits=FALSE, chrommin=NA, chrommax=NA, querymin=NA, querymax=NA) {
  ctginfo <- contiginfo(aligns)
  orderedctginfo <- orderedcontigs(aligns)

  chrom <- aligns[1, "chrom"]
  if(is.na(chrommin)) {
    chrommin <- min(aligns$start)
  }
  if(is.na(chrommax)) {
    chrommax <- max(aligns$end)
  }

  bigcontig <- orderedctginfo[index, "contigname"]
  bigcontigaligns <- aligns[aligns$query==bigcontig, ]

  if (is.na(querymin)) {
    querymin <- ctginfo[ctginfo$contigname==bigcontig, "minpos"]
  }  
  if (is.na(querymax)) {
    querymax <- ctginfo[ctginfo$contigname==bigcontig, "maxpos"]
  }  
  
  contigclusters <- unique(bigcontigaligns$cluster)
  nonsmallcontigclusters <- contigclusters[!str_detect(contigclusters, "Small")]
  
  clusterpalette <- rainbow(length(nonsmallcontigclusters))
  if (length(clusterpalette) == 2) {
    clusterpalette <- c("red", "blue")
  }

  alignclustercolors <- sapply(bigcontigaligns$cluster, function(x) {
    ifelse(str_detect(x, "Small"), "black", clusterpalette[which(nonsmallcontigclusters==x, arr.ind=TRUE)] )
    })
  orderedclusters <- contigclusters[order(contigclusters)]
  legendclustercolors <- sapply(orderedclusters, function(x) {
    ifelse(str_detect(x, "Small"), "black", clusterpalette[which(nonsmallcontigclusters==x, arr.ind=TRUE)] )
    })
  
  plottitle <- ifelse(suppresstitle, "", paste("Alignments of ", genomename, " to benchmark ", benchname))
  xaxtval <- ifelse(suppressxunits, "n", "s")
  xlabval <- ifelse(suppressxunits, "", chrom)
  plot(list(), list(), main=plottitle, xaxt=xaxtval, xlab=xlabval, ylab=orderedctginfo[index, "contigname"], 
       xlim=c(chrommin, chrommax), ylim=c(querymin, querymax))
  segments(x0=bigcontigaligns$start, y0=bigcontigaligns$querystart, x1=bigcontigaligns$end, y1=bigcontigaligns$queryend, col=alignclustercolors)
  legend("bottomright", orderedclusters, pch=15, col=legendclustercolors)
}


multiplotaligns <- function(aligns, chromlength=NA, chromplotfile=NA) {
  querydf <- data.frame(query=unique(aligns[!str_detect(aligns$cluster, "Small"), "query"]))
  querydf$basescovered <-sapply(querydf$query, function(x) {querycoverage(aligns, x)})

  coveredbases <- sum(querydf$basescovered)
  halfchromlength <- 0.5*chromlength/1000000

  if (coveredbases < halfchromlength) {
    return(0)
  }

  if (!is.na(chromplotfile)) {
    pdf(chromplotfile, 8.5, 11.0)
  }

  #par(pardefault)
  numplots <- min(length(querydf$query), 4)

  if (numplots==1) {
    par(mfrow = c(1, 1))
    plotaligns(aligns, 1)
  }
  else {
    orderedquerydf <- querydf[order(querydf$basescovered, decreasing=TRUE), ]
    largestsize <- orderedquerydf[1, "basescovered"]
    smallestsize <- orderedquerydf[length(orderedquerydf$basescovered), "basescovered"]
    
    sizefactor <- largestsize/4.0
    
    plotheights <- ifelse(as.integer(orderedquerydf$basescovered/sizefactor)>=1, as.integer(orderedquerydf$basescovered/sizefactor), 1)
    
    par(mfrow=c(numplots, 1))
    par(oma=c(4,1,3,1))
    nf <- layout(matrix(seq(1, numplots),ncol=1), widths=rep(5, numplots), heights=plotheights, TRUE)
    
    par(mar=c(0,4,2,1))
    plotaligns(aligns, 1, suppressaxis=TRUE)
    if (numplots >= 3) {
       for (plotno in seq(2, numplots-1)) {
         par(mar=c(0,4,0,1))
         plotaligns(aligns, plotno, suppressaxis=TRUE, suppresstitle=TRUE)
       }
    }
    par(mar=c(5,4,0,1))
    plotaligns(aligns, numplots, suppressaxis=FALSE, suppresstitle=TRUE)
  }
  if (!is.na(chromplotfile)) {
    dev.off()
  }

}

querycoverage <- function(aligns, queryname) {
  queryaligns <- aligns[aligns$query==queryname, ]
  covgsum <- sum(abs(queryaligns$queryend - queryaligns$querystart + 1))
  
  return(covgsum)
}

aligns <- readaligns(chromfile)
chromplotfile <- chromfile
chromplotfile <- sub(".bed", ".pdf", chromplotfile)
multiplotaligns(aligns, chromlength=chromlength, chromplotfile)
