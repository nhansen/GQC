import os
import pybedtools
import logging
import math
import pysam
import random
import numpy as np
import array
from GQC import bedtoolslib
from GQC import seqparse

logger = logging.getLogger(__name__)

# Routine choose_bin_size attempts to set the size of bins that will have 
# roughly 1000 primary, non-secondary alignment starts within each bin.
# Input: pysam.AlignmentFile object, pysam.FastaFile object, arguments
# Output: Integer giving the recommended number of bins

def choose_bin_size(alignobj, refobj, includedintervals, args)->int:

    [intervalstart, intervalend] = [0, 0]
    totalreads = 0
    totalrefbases = 0
    for interval in includedintervals:
        ref = interval.chrom
        if ref == "chrM":
            continue
        intervalstart = interval.start
        intervalend = interval.end
        def startsinbin(read):
            if read.reference_start >= intervalstart and read.reference_start < intervalend:
                return True
            else:
                return False

        if args.bincovoverlap:
            readsmapped = alignobj.count(contig = ref, start = intervalstart, end = intervalend)
        else:
            readsmapped = alignobj.count(contig = ref, start = intervalstart, end = intervalend, read_callback = startsinbin)
        totalreads += readsmapped
        totalrefbases += intervalend - intervalstart
        logger.debug("Tallied interval " + ref + ":" + str(intervalstart) + "-" + str(intervalend) + "\t" + str(totalreads) + "\t" + str(totalrefbases))

    arrivalsperbase = totalreads/totalrefbases
    requiredbin = 1000/arrivalsperbase

    # round binsize to nearest 1000 base pairs
    binsize = 1000*math.floor(requiredbin/1000 + 0.5)

    logger.debug("Chose binsize " + str(binsize))

    return binsize
   
def initiate_included_bins(refobj, includedbins, binsize, args)->list:
    bindict = {}
    binnumber = 0

    logger.debug("Initiating bins")
    totallength = 0
    for ref in refobj.references:
        if ref == "chrM":
            continue
        reflength = refobj.get_reference_length(ref)
        totallength += reflength

    numbins = int(totallength/binsize)
    bincounts = array.array('i', (0,)) * numbins
    logger.debug("Zeroed " + str(numbins) + " bins")
    binintervals = []

    #for seqbin in includedbins:
        #if seqbin.end - seqbin.start != binsize:
            #logger.debug("Included interval " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is not complete--skipping")
            #continue
        #else:
            #logger.debug("Included interval " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is complete--writing count")
            #binname = seqbin.chrom + ":" + str(seqbin.start) 
            #binindex = binindexdict[binname]
            #bincount = coveragebincounts[binindex]
            #cfh.write(seqbin.chrom + "\t" + str(seqbin.start) + "\t" + str(seqbin.end) + "\t" + str(bincount) + "\n")

    for ref in refobj.references:
        logger.debug("Binning " + ref)
        if ref == "chrM":
            continue
        reflength = refobj.get_reference_length(ref)
        binstart = 0
        while binstart + binsize <= reflength:
            binintervals.append(ref + "\t" + str(binstart) + "\t" + str(binstart + binsize) + "\n")
            #binintervalstring = binintervalstring + ref + "\t" + str(binstart) + "\t" + str(binstart + binsize) + "\n"
            bindict[ref + ":" + str(binstart)] = binnumber
            binnumber += 1
            binstart += binsize

    logger.debug("Populated index with " + str(binnumber) + " bin names")

    numbins = len(bincounts)
    bedintervalstring = ''.join(binintervals)
    binintervals = pybedtools.BedTool(binintervalstring, from_string = True)
    numints = len(binintervals)

    logger.debug("Zeroed " + str(numbins) + " bins and created a BedTool object with " + str(numints) + " intervals")

    return [binintervals, bincounts, bindict]

def initiate_bins(refobj, binsize, args)->list:

    bindict = {}
    binnumber = 0

    logger.debug("Initiating bins")
    totallength = 0
    for ref in refobj.references:
        if ref == "chrM":
            continue
        reflength = refobj.get_reference_length(ref)
        totallength += reflength

    numbins = int(totallength/binsize)
    bincounts = array.array('i', (0,)) * numbins
    logger.debug("Zeroed " + str(numbins) + " bins")
    binintervals = []

    for ref in refobj.references:
        logger.debug("Binning " + ref)
        if ref == "chrM":
            continue
        reflength = refobj.get_reference_length(ref)
        binstart = 0
        while binstart + binsize <= reflength:
            binintervals.append(ref + "\t" + str(binstart) + "\t" + str(binstart + binsize) + "\n")
            bindict[ref + ":" + str(binstart)] = binnumber
            binnumber += 1
            binstart += binsize

    logger.debug("Populated index with " + str(binnumber) + " bin names")

    numbins = len(bincounts)
    binintervalstring = ''.join(binintervals)
    binintervals = pybedtools.BedTool(binintervalstring, from_string = True)
    numints = len(binintervals)

    logger.debug("Zeroed " + str(numbins) + " bins and created a BedTool object with " + str(numints) + " intervals")

    return [binintervals, bincounts, bindict]

def write_nonexcluded_bin_counts(outputbed:str, includedintervals:pybedtools.BedTool, coveragebins:pybedtools.BedTool, coveragebincounts:list, binindexdict:dict, binsize:int):
    # this should truncate intervals that aren't fully included in the includedintervals:
    includedbins = bedtoolslib.intersectintervals(coveragebins, includedintervals)
    with open(outputbed, "w") as cfh:
        for seqbin in includedbins:
            if seqbin.end - seqbin.start != binsize:
                logger.debug("Included interval " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is not complete--skipping")
                continue
            else:
                logger.debug("Included interval " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is complete--writing count")
                binname = seqbin.chrom + ":" + str(seqbin.start) 
                binindex = binindexdict[binname]
                bincount = coveragebincounts[binindex]
                cfh.write(seqbin.chrom + "\t" + str(seqbin.start) + "\t" + str(seqbin.end) + "\t" + str(bincount) + "\n")

    return 0

def bin_gccontent_extreme_kmers(refobj:pysam.FastaFile, chrom:str, start:int, end:int, kmersize:int):
    logger.debug("Calculating extreme kmer counts for " + chrom + ":" + str(start) + "-" + str(end))
    refsequence = refobj.fetch(reference=chrom, start=start, end=end).upper()
    [count_at, count_ag, count_ac, count_gc] = [0, 0, 0, 0]
    # also track the GC content in each bin
    [count_gcnucs, count_atnucs] = [0, 0]
    seqstart = 0
    seqstop = kmersize
    # count bases in initial kmer:
    basecounts = {'A':0, 'C':0, 'G':0, 'T':0, 'N':0}
    binsize = end - start
    if seqstop < binsize:
        for i in (range(0, seqstop)):
            thisnuc = refsequence[i]
            basecounts[thisnuc] += 1
            if thisnuc == 'G' or thisnuc == 'C':
                count_gcnucs += 1
            else:
                count_atnucs += 1
        kmerseq = refsequence[0:seqstop]

    while seqstop < binsize:
        if basecounts['N']==0:
            if basecounts['G']==0 and basecounts['C']==0:
                count_at += 1
            if basecounts['A']==0 and basecounts['T']==0:
                count_gc += 1
            if (basecounts['T']==0 and basecounts['C']==0) or (basecounts['A']==0 and basecounts['G']==0):
                count_ag += 1
            if (basecounts['T']==0 and basecounts['G']==0) or (basecounts['A']==0 and basecounts['C']==0):
                count_ac += 1

        # subtract the first base and add the next:
        nextbase = refsequence[seqstop]
        basecounts[kmerseq[0]] -= 1
        basecounts[nextbase] += 1
        kmerseq = kmerseq[1:] + nextbase
        if nextbase == 'G' or nextbase == 'C':
            count_gcnucs += 1
        else:
            count_atnucs += 1

        seqstop += 1

    return [count_at, count_ag, count_ac, count_gc, count_gcnucs, count_atnucs]


def count_extreme_kmers_in_bins(outputbed:str, refobj:pysam.FastaFile, includedintervals:pybedtools.BedTool, coveragebins:pybedtools.BedTool, binsize:int):
    # this should truncate intervals that aren't fully included in the includedintervals:
    includedbins = bedtoolslib.intersectintervals(coveragebins, includedintervals)
    kmersize = 40

    with open(outputbed, "w") as cfh:
        for seqbin in includedbins:
            if seqbin.end - seqbin.start != binsize:
                logger.debug("Bin containing " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is partially excluded--skipping")
                continue
            else:
                [count_at, count_ag, count_ac, count_gc, count_gcnucs, count_atnucs] =  bin_gccontent_extreme_kmers(refobj, seqbin.chrom, seqbin.start, seqbin.end, kmersize)
                cfh.write(seqbin.chrom + "\t" + str(seqbin.start) + "\t" + str(seqbin.end) + "\t" + str(count_at) + "\t" + str(count_ag) + "\t" + str(count_ac) + "\t" + str(count_gc) + "\t" + str(count_gcnucs) + "\t" + str(count_atnucs) + "\n")
                logger.debug(seqbin.chrom + "\t" + str(seqbin.start) + "\t" + str(seqbin.end) + "\t" + str(count_at) + "\t" + str(count_ag) + "\t" + str(count_ac) + "\t" + str(count_gc) + "\t" + str(count_gcnucs) + "\t" + str(count_atnucs) + "\n")

    return 0

def tally_chrom_starts_ends(alignobj:pysam.AlignmentFile, refobj:pysam.FastaFile, chrom:str, outputcountbed:str, args):

    chromlength = refobj.get_reference_length(chrom)
    startendcountlist = np.zeros(chromlength, dtype=int)

    for align in alignobj.fetch(chrom):

        refstart = align.reference_start
        refend = align.reference_end

        if align.is_secondary or refstart is None or refend is None:
            continue

        startendcountlist[refstart] += 1
        startendcountlist[refend - 1] -= 1

    print("Done with " + chrom)
    logger.debug("Done with " + chrom)

    return(startendcountlist)

def tally_bin_coverages(alignobj:pysam.AlignmentFile, refobj:pysam.FastaFile, includedintervals:pybedtools.BedTool, outputcountbed:str, outputincludedcountbed, args):

    if args.covbinsize == 0:
        binsize = choose_bin_size(alignobj, refobj, includedintervals, args)
    else:
        binsize = args.covbinsize

    logger.info("Using bin size " + str(binsize))

    # tally coverage across chromosomes, mosdepth style:
    with open(outputcountbed, "w") as cfh:
        for chrom in alignobj.references:
            startendarray = tally_chrom_starts_ends(alignobj, refobj, chrom, outputcountbed, args)
            # binstart/binend is the first/last zero based position in the current bin
            binstart = 0
            binend = binstart + binsize - 1
            currentreads = 0
            # 0-based position
            currentpos = 0
            chromlength = refobj.get_reference_length(chrom)
            bintotal = 0
    
            while currentpos < chromlength:
                # calculate how many reads are covering this position:
                currentreads = currentreads + startendarray[currentpos]
                bintotal += currentreads
                if currentpos == binend:
                    binaverage = 1.0*bintotal/binsize
                    [count_at, count_ag, count_ac, count_gc, count_gcnucs, count_atnucs] =  bin_gccontent_extreme_kmers(refobj, chrom, binstart, binend + 1, 40)
                    cfh.write(chrom + "\t" + str(binstart) + "\t" + str(binend + 1) + "\t" + str(count_gcnucs) + "\t" + str(binaverage) + "\n")
                    bintotal = 0
                    binstart = binend + 1
                    binend = binstart + binsize - 1
                currentpos += 1

            print("Current reads at end of " + chrom + " is " + str(currentreads))
            logger.debug("Current reads at end of " + chrom + " is " + str(currentreads))

    includedcountfile = bedtoolslib.intersectbed(outputcountbed, includedbed, requirewhole=True, writefirst=True)
    logger.info("Wrote count bed file of included whole bins")

def tally_included_bin_coverages(alignobj:pysam.AlignmentFile, refobj:pysam.FastaFile, includedintervals:pybedtools.BedTool, outputcountbed:str, args):

    if args.covbinsize == 0:
        binsize = choose_bin_size(alignobj, refobj, includedintervals, args)
    else:
        binsize = args.covbinsize

    logger.info("Using bin size " + str(binsize))
    [coveragebins, coveragebincounts, coveragebinindex] = initiate_bins(alignobj, binsize, args)

    # this should truncate intervals that aren't fully included in the includedintervals:
    ibin = 0
    kmersize = 40
    includedbins = bedtoolslib.intersectintervals(coveragebins, includedintervals)
    
    [binstart, binend] = [0, 0]
    def startsinbin(read):
        if read.reference_start >= binstart and read.reference_start < binend:
            return True
        else:
            return False

    with open(outputcountbed, "w") as cfh:
        for seqbin in includedbins:
            if seqbin.chrom == "chrM":
                continue
            if seqbin.end - seqbin.start != binsize:
                logger.debug("Included interval " + seqbin.chrom + ":" + str(seqbin.start) + "-" + str(seqbin.end) + " is not complete--skipping")
                continue
            else:
                binchrom = seqbin.chrom
                binstart = seqbin.start
                binend = seqbin.end
                if args.bincovoverlap:
                    bincount = alignobj.count(contig=binchrom, start=binstart, stop=binend)
                else:
                    bincount = alignobj.count(contig=binchrom, start=binstart, stop=binend, read_callback=startsinbin)
                coveragebincounts[ibin] = bincount
                # calculate things about this bin:
                [count_at, count_ag, count_ac, count_gc, count_gcnucs, count_atnucs] =  bin_gccontent_extreme_kmers(refobj, binchrom, binstart, binend, kmersize)
                logger.debug(binchrom + "\t" + str(binstart) + "\t" + str(binend) + "\t" + str(bincount))
                cfh.write(binchrom + "\t" + str(binstart) + "\t" + str(binend) + "\t" + str(count_at) + "\t" + str(count_ag) + "\t" + str(count_ac) + "\t" + str(count_gc) + "\t" + str(count_gcnucs) + "\t" + str(count_atnucs) + "\t" +  str(bincount) + "\n")

def read_extreme_kmer_counts(countfile:str)->dict:

    extremecounts = {}
    with open(countfile, "r") as cfh:
        countline = cfh.readline()
        while countline:
            countline = countline.rstrip()
            [kmertype, count] = countline.split("\t")
            extremecounts[kmertype] = int(count)
            countline = cfh.readline()

    return(extremecounts)

def compare_read_kmers_to_benchmark_kmers(alignobj:pysam.AlignmentFile, refobj:pysam.FastaFile, includedintervals:pybedtools.BedTool, outputbenchmarkcountbed:str, outputreadseqcountbed, outputstrandseqcountbed, config, args):
 
    # Should have this read in from an included count file specified in the config, but for now, calculating it:
    kmersize = args.covkmersize

    if "extremekmercounts" in config.keys() and os.path.exists(config["extremekmercounts"]):
        benchmark_kmer_counts = read_extreme_kmer_counts(config["extremekmercounts"])
    else:
        benchmark_kmer_counts = tally_included_benchmark_kmers(refobj, includedintervals, args)

        with open(outputbenchmarkcountbed, "w") as ofh:
            for kmerseq in benchmark_kmer_counts.keys():
                kmercount = benchmark_kmer_counts[kmerseq]
                ofh.write(kmerseq + "\t" + str(kmercount) + "\n")

    [aligned_read_kmer_counts, stranded_read_kmer_counts] = tally_included_aligned_read_kmers(alignobj, refobj, includedintervals, args)

    with open(outputreadseqcountbed, "w") as ofh:
        for kmerseq in aligned_read_kmer_counts.keys():
            kmercount = aligned_read_kmer_counts[kmerseq]
            if kmerseq in benchmark_kmer_counts.keys():
                benchkmercount = benchmark_kmer_counts[kmerseq]
                kmerratio = 1.0*kmercount/benchkmercount
            else:
                benchkmercount = 'NA'
                kmerratio = 'NA'

            ofh.write(kmerseq + "\t" + str(kmercount) + "\t" + str(benchkmercount) + "\t" + str(kmerratio) + "\n")

    with open(outputstrandseqcountbed, "w") as ofh:
        for kmerseq in stranded_read_kmer_counts.keys():
            kmercount = stranded_read_kmer_counts[kmerseq]
            if kmerseq == 'AT' or kmerseq == 'GC' or kmerseq == 'NA': # own reverse complement, so use total benchmark count
                benchkmercount = benchmark_kmer_counts[kmerseq]
                kmerratio = 1.0*kmercount/benchkmercount
            elif kmerseq in benchmark_kmer_counts.keys():
                benchkmercount = 0.5 * benchmark_kmer_counts[kmerseq]
                kmerratio = 1.0*kmercount/benchkmercount
            else:
                rckmerseq = seqparse.revcomp(kmerseq)
                if rckmerseq in benchmark_kmer_counts.keys():
                    benchkmercount = 0.5 * benchmark_kmer_counts[rckmerseq]
                    kmerratio = 1.0*kmercount/benchkmercount
                else:
                    benchkmercount = 'NA'
                    kmerratio = 'NA'
            ofh.write(kmerseq + "\t" + str(kmercount) + "\t" + str(benchkmercount) + "\t" + str(kmerratio) + "\n")

#
# This function traverses all included intervals (the whole genome or regions specified with --regions, subtracting
# those not in the excluded_regions file), and tallies kmers in the benchmark consensus sequence. chrM is not 
# included due to its unusual copy number.
#
def tally_included_benchmark_kmers(refobj, includedintervals, args)->dict:
    kmercountdict = {}

    k = args.covkmersize
    for interval in includedintervals:
        if interval.chrom == "chrM":
            logger.debug("Skipping mitochondrial region from " + str(interval.start) + " to " + str(interval.end))
            continue
        logger.debug("Calculating kmer counts for " + interval.chrom + ":" + str(interval.start) + "-" + str(interval.end))
        refsequence = refobj.fetch(reference=interval.chrom, start=interval.start, end=interval.end).upper()
        seqstart = 0
        seqstop = k
        # count initial kmer:
        kmerseq = refsequence[0:seqstop]
        # compare to reverse comp and tally
        rckmerseq = seqparse.revcomp(kmerseq)

        intervallength = interval.end - interval.start
        tallykmer = "NA"
        while seqstop <= intervallength:
            if k >= 5: # record extreme kmer types (only 6 different ones since AA==TT, CC==GG, AC==GT and AG==CT for the double stranded benchmark)
                acount = kmerseq.count("A")
                tcount = kmerseq.count("T")
                gcount = kmerseq.count("G")
                ccount = kmerseq.count("C")
                if acount==0 and tcount==0: # kmer is all G's and/or C's
                    if gcount == 0 or ccount == 0: # kmer is entirely G or entirely C
                        tallykmer = "CC"
                    else: # kmer is mix of G's and C's
                        tallykmer = "GC"
                elif gcount==0 and ccount==0: # kmer is all A's and/or T's
                    if acount == 0 or tcount == 0: # kmer is entirely A's or entirely T's
                        tallykmer = "AA"
                    else: # kmer is mix of A's and T's
                        tallykmer = "AT"
                elif (gcount==0 and tcount==0) or (ccount==0 and acount==0): # kmer must be all A's and C's or all G's and T's
                    tallykmer = "AC"
                elif (gcount==0 and acount==0) or (ccount==0 and tcount==0): # kmer must be all A's and G's or all T's and C's
                    tallykmer = "AG"
                else: # kmer has at least three different nucleotides
                    tallykmer = "NA"
            elif kmerseq <= rckmerseq: # if using small kmers (right now, just 3), tally the lower-sorted ("canonical") kmer
                tallykmer = kmerseq
            else:
                tallykmer = rckmerseq
            if tallykmer in kmercountdict.keys():
                kmercountdict[tallykmer] += 1
            else:
                kmercountdict[tallykmer] = 1

            # subtract the first base and add the next:
            if seqstop == intervallength:
                break
            nextbase = refsequence[seqstop]
            kmerseq = kmerseq[1:] + nextbase
            rckmerseq = seqparse.revcomp(kmerseq)
            seqstop += 1

        logger.debug("Last kmer recorded: " + tallykmer)

    return kmercountdict

def tally_included_aligned_read_kmers(alignobj, refobj, includedintervals, args)->dict:
    canonicalkmercountdict = {}
    strandedkmercountdict = {}

    k = args.covkmersize
    minreadalignedpercentage = args.minreadalignedpercentage
    for interval in includedintervals:
        if interval.chrom == "chrM":
            logger.debug("Skipping read alignments in mitochondrial region from " + str(interval.start) + " to " + str(interval.end))
            continue
        logger.debug("Calculating kmer counts for read sequences aligned across " + interval.chrom + ":" + str(interval.start) + "-" + str(interval.end))
        regionstring = interval.chrom + ":" + str(interval.start + 1) + "-" + str(interval.end)

        for align in alignobj.fetch(region=regionstring):
            if align.is_secondary or align.cigartuples is None or (args.downsample is not None and random.random() >= args.downsample):
                continue

            # read sequence is align.query_sequence[query_alignment_start:query_alignment_end], which doesn't include clipped ends
            # it needs to be trimmed to exclude reference positions that are outside of the included region
            alignedreadseq = align.query_alignment_sequence
            querylength = len(alignedreadseq)
            readlength = align.infer_query_length()
            alignrefstart = align.reference_start
            alignrefend = align.reference_end
            if 100*querylength/readlength < minreadalignedpercentage:
                percentclipped = int(100*(readlength - querylength)/readlength)
                logger.debug("Skipping read in " + regionstring + "(" + str(alignrefstart) + "-" + str(alignrefend) + ") due to excessive clipping: " + str(percentclipped) + " percent")
                continue
            cigarops = align.cigartuples
            cigarindex = 0
            while interval.start > alignrefstart:
                logger.debug(str(interval.start) + " is greater than " + str(alignrefstart))
                cigarop = cigarops[cigarindex]
                if cigarop[0] in [0, 1]:
                    alignedreadseq = alignedreadseq[cigarop[1]:]
                if cigarop[0] in [0, 2]:
                    alignrefstart += cigarop[1]
                cigarindex += 1
                if cigarindex >= len(cigarops) - 1:
                    logger.debug("Reached end of cigar ops while trimming alignment start!")
                    break
            cigarindex = len(cigarops) - 1
            while interval.end < alignrefend:
                logger.debug(str(interval.end) + " is less than " + str(alignrefend))
                cigarop = cigarops[cigarindex]
                if cigarop[0] in [0, 1]:
                    alignedreadseq = alignedreadseq[0:-1*cigarop[1]]
                if cigarop[0] in [0, 2]:
                    alignrefend -= cigarop[1]
                cigarindex -= 1
                if cigarindex < 0:
                    logger.debug("Reached beginning of cigar ops while trimming alignment end!")
                    break

            logger.debug("Calculating kmer counts for read " + align.query_name + " from " + align.reference_name + ":" + str(alignrefstart) + "-" + str(alignrefend))
            seqstart = 0
            seqstop = k
            # count initial kmer:
            querylength = len(alignedreadseq)
            kmerseq = alignedreadseq[0:seqstop]
            # compare to reverse comp and tally
            rckmerseq = seqparse.revcomp(kmerseq)
    
            queryreversestrand = align.is_reverse
            tallykmer = "NA"
            while seqstop <= querylength:
                reckmerseq = kmerseq
                recrckmerseq = rckmerseq
                if k >= 5:
                    acount = kmerseq.count("A")
                    tcount = kmerseq.count("T")
                    gcount = kmerseq.count("G")
                    ccount = kmerseq.count("C")
                    if acount==0 and tcount==0: # kmer is all G's and/or C's
                        if gcount==0: # all C's
                            tallykmer = "CC"
                            reckmerseq = "CC"
                            recrckmerseq = "GG"
                        elif ccount == 0: # all G's
                            tallykmer = "CC"
                            reckmerseq = "GG"
                            recrckmerseq = "CC"
                        else: # all G's and C's
                            tallykmer = "GC"
                            reckmerseq = "GC"
                            recrckmerseq = "GC"
                    elif gcount==0 and ccount==0: # kmer is all A's and/or T's
                        if tcount==0: # all A's
                            tallykmer = "AA"
                            reckmerseq = "AA"
                            recrckmerseq = "TT"
                        if acount==0: # all T's
                            tallykmer = "AA"
                            reckmerseq = "TT"
                            recrckmerseq = "AA"
                        else: # all A's and T's
                            tallykmer = "AT"
                            reckmerseq = "AT"
                            recrckmerseq = "AT"
                    elif gcount==0 and tcount==0: # kmer is all A's and C's
                        tallykmer = "AC"
                        reckmerseq = "AC"
                        recrckmerseq = "GT"
                    elif ccount==0 and acount==0: # kmer is all G's and T's
                        tallykmer = "AC"
                        reckmerseq = "GT"
                        recrckmerseq = "AC"
                    elif ccount==0 and tcount==0: # kmer is all A's and G's
                        tallykmer = "AG"
                        reckmerseq = "AG"
                        recrckmerseq = "CT"
                    elif gcount==0 and acount==0: # kmer is all C's and T's
                        tallykmer = "AG"
                        reckmerseq = "CT"
                        recrckmerseq = "AG"
                    else:
                        tallykmer = "NA"
                        reckmerseq = "NA"
                        recrckmerseq = "NA"
                elif kmerseq <= rckmerseq:
                    tallykmer = kmerseq
                else:
                    tallykmer = rckmerseq

                # record the canonical (strand independent) kmer
                if tallykmer in canonicalkmercountdict.keys():
                    canonicalkmercountdict[tallykmer] += 1
                else:
                    canonicalkmercountdict[tallykmer] = 1

                if queryreversestrand:
                    strandedkmer = recrckmerseq
                else:
                    strandedkmer = reckmerseq
                if strandedkmer in strandedkmercountdict.keys():
                    strandedkmercountdict[strandedkmer] += 1
                else:
                    strandedkmercountdict[strandedkmer] = 1
    
                # subtract the first base and add the next:
                if seqstop == querylength:
                    break
                nextbase = alignedreadseq[seqstop]
                kmerseq = kmerseq[1:] + nextbase
                rckmerseq = seqparse.revcomp(kmerseq)
                seqstop += 1

            logger.debug("Last kmer recorded: " + tallykmer + "/" + kmerseq)

    return [canonicalkmercountdict, strandedkmercountdict]

