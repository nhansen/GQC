library(stringr)

#setwd("/Users/nhansen/OneDrive/HPRC_assembly_comparison/data_for_Rplots/alignplots/hprc_minimap2")
args = commandArgs(trailingOnly=TRUE)

pardefault <- par()

chromname <- ifelse(!( is.na(args[1])), args[1], "HG01346#1#CM086596.1")
genomename <- ifelse(!( is.na(args[2])), args[2], "HG01346_HPRC")
benchname <- ifelse(!( is.na(args[3])), args[3], "HG01346_verkko_curated_y")
outputdir <- ifelse(!( is.na(args[4])), args[4], ".")
chromlength <- ifelse(! ( is.na(args[5])), as.integer(args[5]), 50304960)
refgapfile <- ifelse(! ( is.na(args[6])), args[6], NA)
querygapfile <- ifelse(! ( is.na(args[7])), args[7], NA)
vertline <- ifelse(! ( is.na(args[8])), as.integer(args[8]), NA)

# plot positions in megabases:
axisunitbases <- 1000000
safe_colorblind_palette <- c("#332288", "#88CCEE", "#CC6677", "#DDCC77", "#117733", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888",
			     "#648fff", "#fe6100", "#785ef0")

readaligns <- function(chromfile) {
  aligns <- read.table(chromfile, sep="\t", header=FALSE, comment.char="!")
  names(aligns) <- c("chrom", "start", "end", "alignname")

  aligns$start <- as.numeric(aligns$start)/axisunitbases
  aligns$end <- as.numeric(aligns$end)/axisunitbases
  aligns$query <- sapply(seq(1, length(aligns$alignname)), function(i) {us_fields <- strsplit(aligns$alignname, split="_"); return(paste(us_fields[[i]][1:(length(us_fields[[i]])-3)], sep="_", collapse="_"))})
  aligns$querystart <- sapply(seq(1, length(aligns$alignname)), function(i) {us_fields <- strsplit(aligns$alignname, split="_"); return(as.numeric(us_fields[[i]][length(us_fields[[i]]) - 2])/axisunitbases)})
  aligns$queryend <- sapply(seq(1, length(aligns$alignname)), function(i) {us_fields <- strsplit(aligns$alignname, split="_"); return(as.numeric(us_fields[[i]][length(us_fields[[i]]) - 1])/axisunitbases)})
  aligns$cluster <- sapply(seq(1, length(aligns$alignname)), function(i) {us_fields <- strsplit(aligns$alignname, split="_"); return(us_fields[[i]][length(us_fields[[i]])])})

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
  orderedctginfo <- ctginfo[order(ctginfo$minrefpos, decreasing = TRUE), ]
  
  return(orderedctginfo)
}


multiplotaligns <- function(aligns, chromlength=NA, chromplotfile=NA, minfrac=0.5, refgapfile=NA, querygapfile=NA, vertlinepos=NA, maxqueryentries=10, suppresstitle=FALSE, chrommin=NA, chrommax=NA, querymin=NA, querymax=NA, unitbases=axisunitbases, querytickdistance=10) {

  chromname <- unique(aligns[1]["chrom"])  
  # largeclusterqueryaligninfo dataframe has only the contig names for contigs in the large clusters of alignments:
  largeclusterqueryaligninfo <- data.frame(query=unique(aligns[!str_detect(aligns$cluster, "Small"), "query"]))
  # querybasescovered is the total aligned bases in each contig (including those aligned within small clusters)
  largeclusterqueryaligninfo$querybasescovered <-sapply(largeclusterqueryaligninfo$query, function(x) {querycoverage(aligns, x)})
  largeclusterqueryaligninfo$refbasescovered <-sapply(largeclusterqueryaligninfo$query, function(x) {refcoverage(aligns, x)})
  largeclusterqueryaligninfo$minrefpos <-sapply(largeclusterqueryaligninfo$query, function(x) {minrefpos(aligns, x)})
  largeclusterqueryaligninfo$minquerypos <-sapply(largeclusterqueryaligninfo$query, function(x) {minpos(aligns, x)})
  largeclusterqueryaligninfo$maxquerypos <-sapply(largeclusterqueryaligninfo$query, function(x) {maxpos(aligns, x)})  # probably want to added some kind of position-weighted average ref coordinate measure
  largeclusterqueryaligninfo$queryaxislength <-sapply(largeclusterqueryaligninfo$query, function(x) {largeclusterqueryaligninfo[largeclusterqueryaligninfo$query==x, "maxquerypos"]-largeclusterqueryaligninfo[largeclusterqueryaligninfo$query==x, "minquerypos"]})

  # don't bother to plot any alignments if total aligned query bases are less than minfrac (0.5) of the ref chrom length:
  coveredbases <- sum(largeclusterqueryaligninfo$querybasescovered)
  halfchromlength <- 0.5*chromlength/unitbases
  if (coveredbases < halfchromlength) {
    return(0)
  }

  # open the file to plot to if there is one:
  if (!is.na(chromplotfile)) {
    pdf(chromplotfile, 8.5, 11.0)
  }

  # read in gaps if there are gap files:
  refgaps <- gaplocs(refgapfile)
  querygaps <- gaplocs(querygapfile)

  # store graph parameters so they can be restored after changing:
  par(pardefault)
  
  # order the contigs along the reference chromosome and calculate pseudocoordinates:
  positionorderedqueryinfo <- largeclusterqueryaligninfo[order(largeclusterqueryaligninfo$minrefpos, decreasing=FALSE), ]
  numsections <- min(length(largeclusterqueryaligninfo$query), maxqueryentries)
  positionorderedqueryinfo$pseudostart <- sapply(seq(1, length(positionorderedqueryinfo$query)), function(x) { if (x==1) {return(0)} else {return(sum(positionorderedqueryinfo[1:(x-1), "queryaxislength"]))} })
  positionorderedqueryinfo$pseudoend <- sapply(seq(1, length(positionorderedqueryinfo$query)), function(x) { return(sum(positionorderedqueryinfo[(1:x), "queryaxislength"])) })

  bigcontiglist <- unique(positionorderedqueryinfo[ , "query"])
  bigcontigaligns <- lapply(positionorderedqueryinfo$query, function(x) {return(aligns[aligns$query==x, ])})

  bigcontigalignlist <- aligns[is.element(aligns$query, bigcontiglist),]
  contigclusters <- unique(bigcontigalignlist$cluster)
  nonsmallcontigclusters <- contigclusters[!str_detect(contigclusters, "Small")]
  
  clusterpalette <- safe_colorblind_palette[1:length(nonsmallcontigclusters)]
  bigclustercolors <- data.frame("clustername"=nonsmallcontigclusters, "clustercolor"=clusterpalette)
 
  orderedclusters <- str_sort(contigclusters, numeric=TRUE)
 
  #orderedclusters <- contigclusters[order(contigclusters)]
  legendclustercolors <- sapply(orderedclusters, function(x) {
    ifelse(str_detect(x, "Small"), "black", clusterpalette[which(nonsmallcontigclusters==x, arr.ind=TRUE)] )
  })
  
  plottitle <- ifelse(suppresstitle, "", paste("Alignments of ", genomename, " to ", benchname))
  xlabval <- paste(c(chromname, " (Mb)"), sep="", collapse="" )
  if (is.na(chrommin)) {
    chrommin <- min(positionorderedqueryinfo$minrefpos)
  }
  if (is.na(chrommax)) {
    chrommax <- chromlength/unitbases
  }
  if (is.na(querymin)) {
    querymin <- 0
  }
  if (is.na(querymax)) {
    querymax <- max(positionorderedqueryinfo$pseudoend + 1)
  }
  
  plot(list(), list(), main=plottitle, xaxs='i', yaxt='n', xlab=xlabval, yaxs='i', ylab="", 
       xlim=c(chrommin, chrommax), ylim=c(querymin, querymax))

  contigdividelines <- positionorderedqueryinfo$pseudostart
  abline(h=contigdividelines[-1])

  alltickpseudopositions = c()
  allticklabels = c()
  for(contigaligns in bigcontigaligns) {
    contigname <- unique(contigaligns$query)[[1]]
    contigmin <- positionorderedqueryinfo[positionorderedqueryinfo$query==contigname, "minquerypos"]
    contigmax <- positionorderedqueryinfo[positionorderedqueryinfo$query==contigname, "maxquerypos"]
    contigpseudostart <- positionorderedqueryinfo[positionorderedqueryinfo$query==contigname, "pseudostart"]
    contigpseudoend <- positionorderedqueryinfo[positionorderedqueryinfo$query==contigname, "pseudoend"]
    if ((ceiling(contigmin/querytickdistance) < contigmax/querytickdistance) && (floor(contigmax/querytickdistance) > contigmin/querytickdistance)) {
      contigtickpositions <- querytickdistance*seq.int(ceiling(contigmin/querytickdistance), floor(contigmax/querytickdistance))
    }
    else {
      contigtickpositions <- c()
    }
    contigtickpseudopositions <- contigtickpositions + contigpseudostart - contigmin
    contigticklabels <- as.character(contigtickpositions)
    alltickpseudopositions <- c(alltickpseudopositions, contigtickpseudopositions)
    allticklabels <- c(allticklabels, contigticklabels)
    segmentx0vals <- contigaligns$start
    segmentx1vals <- contigaligns$end
    segmenty0vals <- sapply(seq(1, length(contigaligns$query)), function(x) { contigaligns[x, "querystart"] + contigpseudostart - contigmin})
    segmenty1vals <- sapply(seq(1, length(contigaligns$query)), function(x) { contigaligns[x, "queryend"] + contigpseudostart - contigmin})
    segmentcolors <- sapply(contigaligns$cluster, function(x) {
      ifelse(str_detect(x, "Small"), "black", bigclustercolors[which(bigclustercolors$clustername==x, arr.ind=TRUE), "clustercolor"] )
    })

    cexfactor <- 0.8 # will want to adjust this with the lengths of the contig names
    contiglabelheight <- (contigpseudostart + contigpseudoend)/2.0
    if ((querymin != 0) & ((contigpseudostart <= querymin) | (contigpseudoend >= querymax))) {
      # adjust height of label
      contiglabelheight <- ((querymin + querymax)/2.0)
    }
    #return(contiglabelheight)
    if (contigpseudoend - contigpseudostart > 10) {
      mtext(contigname, side=2, line=2, at=contiglabelheight, cex=cexfactor)
    }
    
    segments(x0=segmentx0vals, y0=segmenty0vals, x1=segmentx1vals, y1=segmenty1vals, col=segmentcolors, lwd=1.5)
  }
  # ticks/labels:
  
  inboundarytickmarks <- which(alltickpseudopositions>querymin & alltickpseudopositions<querymax)
  axis(side=2, at=alltickpseudopositions, labels = allticklabels)

  if (!(is.na(vertlinepos))) {
    abline(v=vertlinepos, lty=4)
    if (!is.na(vertlinelabel)) {
      text(vertlinepos+2, 1, labels= vertlinelabel)
    }
  }
  if (length(refgaps)>1 || !(is.na(refgaps))) {
    refstarts <- refgaps[refgaps$contig==chrom, "start"]/unitbases
    refends <- refgaps[refgaps$contig==chrom, "end"]/unitbases
    if (length(refstarts) > 0) {
      rect(refstarts, querymin, refends, querymax, col="gray")
    }
  }
  if (length(querygaps)>1 || !(is.na(querygaps))) {
    querystarts <- querygaps[querygaps$contig==chrom, "start"]/unitbases
    queryends <- querygaps[querygaps$contig==chrom, "end"]/unitbases
    if (length(querystarts) > 0) {
      rect(chrommin, querystarts, chrommax, queryends, col="gray")
    }
  }
  
  legend("bottomright", orderedclusters, pch=15, col=legendclustercolors)
  
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

chromfiles <- list.files(outputdir, pattern = paste(c('.', chromname, '.clusters.bed'), sep="", collapse=""))
aligns <- c()
for (file in chromfiles) {
  aligns <- rbind(aligns, readaligns(paste(c(outputdir, "/", file), sep="", collapse="")))
}
chromplotfile <- paste(outputdir, "/", chromname, ".clustered_aligns.pdf", sep="", collapse="")

if (is.na(vertline)) {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile) 
} else {
  multiplotaligns(aligns, chromlength=chromlength, chromplotfile, refgapfile=refgapfile, querygapfile=querygapfile, yqhstart=vertline)  
}

