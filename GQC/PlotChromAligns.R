library(stringr)

#setwd("/Users/nhansen/OneDrive/HPRC_assembly_comparison/data_for_Rplots/alignplots")
args = commandArgs(trailingOnly=TRUE)

pardefault <- par()

chromfile <- ifelse(!( is.na(args[1])), args[1], "clustered_aligns.chr1_MATERNAL.clusters.bed")
genomename <- ifelse(!( is.na(args[2])), args[2], "year1pat")
benchname <- ifelse(!( is.na(args[3])), args[3], "v1.0.1")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")
chromlength <- ifelse(! ( is.na(args[5])), as.integer(args[5]), NA)
refgapfile <- ifelse(! ( is.na(args[6])), args[6], NA)
querygapfile <- ifelse(! ( is.na(args[7])), args[7], NA)
vertline <- ifelse(! ( is.na(args[8])), as.integer(args[8]), NA)

# plot positions in megabases:
axisunitbases <- 1000000
safe_colorblind_palette <- c("#332288", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")

readaligns <- function(chromfile) {
  aligns <- read.table(chromfile, sep="\t", header=FALSE, comment.char="!")
  names(aligns) <- c("chrom", "start", "end", "alignname")

  aligns$start <- as.numeric(aligns$start)/axisunitbases
  aligns$end <- as.numeric(aligns$end)/axisunitbases
  aligns$query <- sapply(seq(1, length(aligns$alignname)), function(i) {strsplit(aligns$alignname, split="_")[[i]][1]})
  aligns$querystart <- sapply(seq(1, length(aligns$alignname)), function(i) {as.numeric(strsplit(aligns$alignname, split="_")[[i]][2])/axisunitbases})
  aligns$queryend <- sapply(seq(1, length(aligns$alignname)), function(i) {as.numeric(strsplit(aligns$alignname, split="_")[[i]][3])/axisunitbases})
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

minrefpos <- function(aligns, contigname) {
  startsends <- c(aligns[aligns$query==contigname, "start"], aligns[aligns$query==contigname, "end"])
  return(min(startsends))
}

maxrefpos <- function(aligns, contigname) {
  startsends <- c(aligns[aligns$query==contigname, "start"], aligns[aligns$query==contigname, "end"])
  return(max(startsends))
}

contiginfo <- function(aligns) {
  
  ctginfo <- data.frame("contigname"=unique(aligns$query))
  ctginfo$minpos <- sapply(ctginfo$contigname, function(x) {minpos(aligns, x)})
  ctginfo$maxpos <- sapply(ctginfo$contigname, function(x) {maxpos(aligns, x)})
  ctginfo$minrefpos <- sapply(ctginfo$contigname, function(x) {minrefpos(aligns, x)})
  ctginfo$maxrefpos <- sapply(ctginfo$contigname, function(x) {maxrefpos(aligns, x)})
  ctginfo$basescovered <- ctginfo$maxpos - ctginfo$minpos
  
  return(ctginfo)
}


orderedcontigs <- function(aligns) {
  ctginfo <- contiginfo(aligns) 
  orderedctginfo <- ctginfo[order(ctginfo$minrefpos, decreasing = FALSE), ]
  
  return(orderedctginfo)
}

plotaligns <- function(aligns, index=1, suppressaxis=FALSE, suppresstitle=FALSE, suppressxunits=FALSE, chrommin=NA, chrommax=NA, querymin=NA, querymax=NA, vertline=NA, vertlinelabel=NA, refgaps=NA, querygaps=NA) {
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
  
  clusterpalette <- safe_colorblind_palette[1:length(nonsmallcontigclusters)]
  
  alignclustercolors <- sapply(bigcontigaligns$cluster, function(x) {
    ifelse(str_detect(x, "Small"), "black", clusterpalette[which(nonsmallcontigclusters==x, arr.ind=TRUE)] )
    })
  orderedclusters <- contigclusters[order(contigclusters)]
  legendclustercolors <- sapply(orderedclusters, function(x) {
    ifelse(str_detect(x, "Small"), "black", clusterpalette[which(nonsmallcontigclusters==x, arr.ind=TRUE)] )
    })
  
  plottitle <- ifelse(suppresstitle, "", paste("Alignments of ", genomename, " to ", benchname))
  xaxtval <- ifelse(suppressxunits, "n", "s")
  xlabval <- ifelse(suppressxunits, "", chrom)
  plot(list(), list(), main=plottitle, xaxs='i', xaxt=xaxtval, xlab=xlabval, yaxs='i', ylab=orderedctginfo[index, "contigname"], 
       xlim=c(chrommin, chrommax), ylim=c(querymin, querymax))
  #plot(list(), list(), main=plottitle, xaxt=xaxtval, xlab=xlabval, ylab=orderedctginfo[index, "contigname"], 
       #xlim=c(chrommin, chrommax), ylim=c(querymin, querymax))
  

  segments(x0=bigcontigaligns$start, y0=bigcontigaligns$querystart, x1=bigcontigaligns$end, y1=bigcontigaligns$queryend, col=alignclustercolors)
  if (!(is.na(vertline))) {
    abline(v=vertline, lty=4)
    if (!is.na(vertlinelabel)) {
      text(vertline+2, 1, labels= vertlinelabel)
    }
  }
  if (length(refgaps$contig)==1 && (is.na(refgaps))) {
    refstarts <- NA
    refends <- NA
  }
  else {
    refstarts <- refgaps[refgaps$contig==chrom, "start"]/axisunitbases
    refends <- refgaps[refgaps$contig==chrom, "end"]/axisunitbases
    if (length(refstarts) > 0) {
      rect(refstarts, querymin, refends, querymax, col="gray")
    }
  }
  if (length(querygaps$contig)==1 && (is.na(querygaps))) {
    querystarts <- NA
    queryends <- NA
  } else {
    querystarts <- querygaps[querygaps$contig==chrom, "start"]/axisunitbases
    queryends <- querygaps[querygaps$contig==chrom, "end"]/axisunitbases
    if (length(querystarts) > 0) {
      rect(chrommin, querystarts, chrommax, queryends, col="gray")
    }
  }

  legend("bottomright", orderedclusters, pch=15, col=legendclustercolors)
}


multiplotaligns <- function(aligns, chromlength=NA, chromplotfile=NA, plotheights=c(), minfrac=0.5, refgapfile=NA, querygapfile=NA, yqhstart=NA) {
  # query dataframe has only the contig names for contigs in the large clusters of alignments:
  querydf <- data.frame(query=unique(aligns[!str_detect(aligns$cluster, "Small"), "query"]))
  # querybasescovered is the total aligned bases in each contig (including those aligned within small clusters)
  querydf$querybasescovered <-sapply(querydf$query, function(x) {querycoverage(aligns, x)})
  querydf$refbasescovered <-sapply(querydf$query, function(x) {refcoverage(aligns, x)})
  querydf$minrefpos <-sapply(querydf$query, function(x) {minrefpos(aligns, x)})
  
  # don't bother to plot any alignments when aligned query bases are less than minfrac (0.5) of the ref chrom length:
  coveredbases <- sum(querydf$querybasescovered)
  halfchromlength <- 0.5*chromlength/axisunitbases

  if (coveredbases < halfchromlength) {
    return(0)
  }

  # open the file to plot to if there is one:
  if (!is.na(chromplotfile)) {
    pdf(chromplotfile, 8.5, 8.5)
  }

  # if yqhstart is set, plot the start of the Yqh region  
  if (!is.na(yqhstart)) {
    yqhstart <- yqhstart/axisunitbases
  }
  
  # read in gaps if there are gap files:
  refgaps <- gaplocs(refgapfile)
  querygaps <- gaplocs(querygapfile)

  par(pardefault)
  numplots <- min(length(querydf$query), 4)

  # if just one contig aligns, no need for a matrix layout:
  if (numplots==1) {
    par(mfrow = c(1, 1))
    plotaligns(aligns, 1, chrommin=0, chrommax=floor(chromlength/axisunitbases+1), vertline=yqhstart, refgaps=refgaps, querygaps=querygaps)
  }
  else {
    sizeorderedquerydf <- querydf[order(querydf$refbasescovered, decreasing=TRUE), ]
    largestsize <- sizeorderedquerydf[1, "refbasescovered"]
    smallestsize <- sizeorderedquerydf[length(sizeorderedquerydf$refbasescovered), "refbasescovered"]
    
    sizefactor <- largestsize/50.0

    orderedquerydf <- querydf[order(querydf$minrefpos, decreasing=FALSE), ]
    if (length(plotheights)==0) {
      plotheights <- ifelse(as.integer(orderedquerydf$refbasescovered/sizefactor)>=1, as.integer(orderedquerydf$refbasescovered/sizefactor), 1)
      plotheights[length(plotheights)] <- plotheights[length(plotheights)] + 20
    }

    par(mfrow=c(numplots, 1))
    par(oma=c(4,1,3,1))
    nf <- layout(matrix(seq(1, numplots),ncol=1), widths=rep(80, numplots), heights=plotheights, TRUE)
    
    par(mar=c(0,4,2,1))
    plotaligns(aligns, 1, suppressaxis=TRUE, chrommin=0, chrommax=floor(chromlength/axisunitbases+1), suppressxunits=TRUE, vertline=yqhstart, refgaps=refgaps, querygaps=querygaps)
    if (numplots >= 3) {
       for (plotno in seq(2, numplots-1)) {
         par(mar=c(0,4,0,1))
         plotaligns(aligns, plotno, suppressaxis=TRUE, suppressxunits=TRUE, suppresstitle=TRUE, chrommin=0, chrommax=floor(chromlength/axisunitbases+1), vertline=yqhstart, refgaps=refgaps, querygaps=querygaps)
       }
    }
    par(mar=c(5,4,0,1))
    plotaligns(aligns, numplots, suppressaxis=FALSE, suppresstitle=TRUE, chrommin=0, chrommax=floor(chromlength/axisunitbases+1), vertline=yqhstart, vertlinelabel='Yqh', refgaps=refgaps, querygaps=querygaps)
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

refcoverage <- function(aligns, queryname) {
  queryaligns <- aligns[aligns$query==queryname, ]
  covgsum <- sum(abs(queryaligns$end - queryaligns$start + 1))
  
  return(covgsum)
}

gaplocs <- function(bedfile) {

  if (!is.na(bedfile)) {
    gaps <- read.table(bedfile, sep="\t", header=FALSE, comment.char="!")
    names(gaps) <- c("contig", "start", "end", "gapname")    
  } else {
    gaps <- NA
  }
  
  return(gaps)
}

aligns <- readaligns(chromfile)
chromplotfile <- chromfile
chromplotfile <- sub(".bed", ".pdf", chromplotfile)

if (is.na(vertline)) {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile) 
} else {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile, yqhstart=vertline)  
}

