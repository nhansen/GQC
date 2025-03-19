import sys
import os
import re
import shutil
import pysam
import argparse
import logging
import pybedtools
import importlib.resources
from pathlib import Path
from GQC import seqparse
from GQC import alignparse
#from GQC import assemblygraph
from GQC import bedtoolslib
from GQC import heatmap

logger = logging.getLogger(__name__)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print all heterozygous variants in a specified set of assembly alignments"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('--bam1', required=False, help='bam file of alignments in one direction from one haplotype mapped to the other')
    parser.add_argument('--bam2', required=False, help='bam file of alignments in the other direction')
    parser.add_argument('--ref1', type=str, required=False, help='(indexed) reference fasta file for the first bam file')
    parser.add_argument('--ref2', type=str, required=False, help='(indexed) reference fasta file for the second bam file')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=10000, help='minimum length of alignment required to be included in het site gathering')
    parser.add_argument('--windowsize', '--hetwindowsize', type=int, required=False, default=100000, help='window size for reporting heterozygosity levels')
    parser.add_argument('-n', '--non1to1', action='store_true', required=False, help='relax 1to1 (default) requirement of matching alignment endpoints')
    parser.add_argument('--heatmap', action='store_true', required=False, help='print a bed file with windowed variant counts and colors representing heterozygosity levels (red high yellow low)')
    parser.add_argument('--heatmapmincount', type=int, required=False, help='window variant count to correspond to the yellow end of the heatmap spectrum')
    parser.add_argument('--heatmapmaxcount', type=int, required=False, help='window variant count to correspond to the red end of the heatmap spectrum')
    parser.add_argument('-g', '--gfa', type=str, required=False, help='gfa file for assembly graph to guide selection of het alignments')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="hetsites", help='prefix to use in output filenames')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

# reads aligns from a bam file, and returns (a) a BED object with all same-autosome (or X/Y) aligned reference intervals (unmerged)
# (b) a BED object with query intervals, and (c) a dictionary of alignments query-able by alignment name (which is query:start-end/ref:start-end
def read_aligns(bamobj, args):
    aligns = {}

    querybedstring = ''
    refbedstring = ''
    #alignnumber = 1
    for align in bamobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= args.minalignlength:
            query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
            querycoords = query + ":" + str(querystart) + "-" + str(queryend)
            refcoords = ref + ":" + str(refstart) + "-" + str(refend)
            querychrom = query.split("_")[0]
            refchrom = ref.split("_")[0]

            if querychrom == refchrom or (querychrom == "chrX" and refchrom == "chrY") or (querychrom == "chrY" and refchrom == "chrX"):
                alignname = querycoords + "/" + refcoords
                aligns[alignname] = align
                querybedstring = querybedstring + query + "\t" + str(querystart) + "\t" + str(queryend) + "\t" + alignname + "\n"
                refbedstring = refbedstring + ref + "\t" + str(refstart) + "\t" + str(refend) + "\t" + alignname + "\n"

    refalignintervals = pybedtools.BedTool(refbedstring, from_string=True)
    queryalignintervals = pybedtools.BedTool(querybedstring, from_string=True)

    return refalignintervals, queryalignintervals, aligns

def find_corresponding_alignment_pairs(alignobj1, alignobj2, refobj1, refobj2, aligndict1, aligndict2, querybed1, querybed2, args):

    alignmapping = {}
    orthologymapstring1 = ""
    orthologymapstring2 = ""
    # look at aligns in the first direction, e.g., maternal aligned to paternal
    for align1 in alignobj1.fetch():
        query1, querystart1, queryend1, ref1, refstart1, refend1, strand1 = alignparse.retrieve_align_data(align1)
        refcoords1 = ref1 + ":" + str(refstart1) + "-" + str(refend1)
        querycoords1 = query1 + ":" + str(querystart1) + "-" + str(queryend1)
        align1name = querycoords1+"/"+refcoords1
        align2name = refcoords1+"/"+querycoords1

        # when there are reciprocal best alignments between haplotypes with matching endpoints in bam1 and bam2, use them:
        if align1name in aligndict1.keys() and align2name in aligndict2.keys():
            align2 = aligndict2[align2name]
            logger.debug("MATCH\t" + align1name + " and " + align2name)
            query2, querystart2, queryend2, ref2, refstart2, refend2, strand2 = alignparse.retrieve_align_data(align2)
            refcoords2 = ref2 + ":" + str(refstart2) + "-" + str(refend2)
            if querycoords1 == refcoords2 and refend1 - refstart1 >= args.minalignlength:
                alignmapping[align1name] = {'align':align1, 'intervals':None}
                alignmapping[align2name] = {'align':align2, 'intervals':None}
                orthologymapstring1 = orthologymapstring1 + ref1 + "\t" + str(refstart1 - 1) + "\t" + str(refend1) + "\t" + align1name + ".1to1" + "\n"
                orthologymapstring2 = orthologymapstring2 + ref2 + "\t" + str(refstart2 - 1) + "\t" + str(refend2) + "\t" + align2name + ".1to1" + "\n"
            elif refend1 - refstart1 >= args.minalignlength:
                logger.warn("Matching alignment " + align2name + " has unexpected reference coords")
        # without 1-1 requirement, will also allow alignments contained in larger opposite-direction alignments to be included
        elif args.non1to1 and align1name in aligndict1.keys():
            refinterval = pybedtools.BedTool(ref1 + "\t" + str(refstart1) + "\t" + str(refend1), from_string = True)
            intersectingaligns = bedtoolslib.intersectintervals(refinterval, querybed2, wb=True)
            for aligninterval in intersectingaligns:
                align2name = aligninterval[6]
                align2chrom = aligninterval.chrom
                align2start = aligninterval.start
                align2end = aligninterval.end
                logger.debug(align1name + " intersectalign " + align2name)
                querycoords2 = aligninterval.chrom + ":" + str(aligninterval.start) + "-" + str(aligninterval.end)
                if align2name in aligndict2.keys():
                    align2 = aligndict2[align2name]
                    query2, querystart2, queryend2, ref2, refstart2, refend2, strand2 = alignparse.retrieve_align_data(align2)
                    querycoords2 = query2 + ":" + str(querystart2) + "-" + str(queryend2)
                    refcoords2 = ref2 + ":" + str(refstart2) + "-" + str(refend2)
                    # trim to common within endpoints from bam1:
                    [commonref1start, commonref1end] = alignparse.compare_alignments(align1, align2, switched=True)
                    logger.debug("MISMATCH\t" + refcoords1 + "/" + querycoords1 + "\t" + str(commonref1start) + "\t" + str(commonref1end))
                    if commonref1start is not None and commonref1end is not None and commonref1start != -1:
                        if align1name not in alignmapping.keys():
                            alignmapping[align1name] = {'align':align1, 'intervals':[[commonref1start, commonref1end]]}
                        elif alignmapping[align1name]['intervals'] is None:
                            logger.debug("Align " + align1name + " already exists in alignmapping dictionary with no restriction intervals")
                            alignmapping[align1name]['intervals'] = [[commonref1start, commonref1end]]
                        else:
                            alignmapping[align1name]['intervals'].append([commonref1start, commonref1end])
                        orthologymapstring1 = orthologymapstring1 + ref1 + "\t" + str(refstart1 - 1) + "\t" + str(refend1) + "\t" + align1name + "." + str(commonref1start) + "_" + str(commonref1end) + "\n"
                    elif commonref1start is None and queryend1 - querystart1 > 1000000:
                        logger.debug("LONG MISMATCH length " + str(queryend1 - querystart1))
                    # now trim to common within endpoints from bam2:
                    [commonref2start, commonref2end] = alignparse.compare_alignments(align2, align1, switched=True)
                    logger.debug("MISMATCH\t" + refcoords2 + "/" + querycoords2 + "\t" + str(commonref2start) + "\t" + str(commonref2end))
                    if commonref2start is not None and commonref2end is not None and commonref2start != -1:
                        if align2name not in alignmapping.keys():
                            alignmapping[align2name] = {'align':align2, 'intervals':[[commonref2start, commonref2end]]}
                        elif alignmapping[align2name]['intervals'] is None:
                            logger.debug("Align " + align2name + " already exists in alignmapping dictionary with no restriction intervals")
                            alignmapping[align1name]['intervals'] = [[commonref1start, commonref1end]]
                        else:
                            alignmapping[align2name]['intervals'].append([commonref2start, commonref2end])
                        orthologymapstring2 = orthologymapstring2 + ref2 + "\t" + str(refstart2 - 1) + "\t" + str(refend2) + "\t" + align2name + "." + str(commonref2start) + "_" + str(commonref2end) + "\n"
                    elif commonref2start is None and queryend2 - querystart2 > 1000000:
                        logger.debug("LONG MISMATCH length " + str(queryend2 - querystart2))

    with open(args.prefix + ".orthoregions1.bed", "w") as ofh1:
        ofh1.write(orthologymapstring1)

    with open(args.prefix + ".orthoregions2.bed", "w") as ofh2:
        ofh2.write(orthologymapstring2)

    return alignmapping

def find_hets_and_coveredregions(bamobj, refobj, queryobj, alignmap:dict, args):

    # strings which will be used to create bed files with pybedtools:
    aligncoveredregions = ""
    aligncoveredwindows = ""
    hetvariants = ""
    # to avoid duplicates due to overlapping alignments:
    includedvariantdict = {}
    allwindowintervalswithcounts = pybedtools.BedTool("", from_string=True)
    for align in bamobj.fetch():
        query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
        querycoords = query + ":" + str(querystart) + "-" + str(queryend)
        refcoords = ref + ":" + str(refstart) + "-" + str(refend)
        alignname = querycoords + "/" + refcoords

        windowintervalswithcounts = pybedtools.BedTool("", from_string=True)
        boundarylist = []
        if alignname in alignmap.keys():
            if alignmap[alignname]['intervals'] is None:
                boundarylist.append([refstart, refend])
            else:
                for intervalpair in alignmap[alignname]['intervals']:
                    boundarylist.append(intervalpair)

            for intervalpair in boundarylist:
                [restrictedstart, restrictedend] = intervalpair
                # add query\tstart\tstop\tname\tscore\tstrand to covered string:
                if restrictedstart < 1 or restrictedend < 1:
                    logger.debug("Ignoring alignment from " + querycoords + " to " + refcoords + " due to invalid restrictedstart " + str(restrictedstart))
                else:
                    coveredstring = ref + "\t" + str(restrictedstart - 1) + "\t" + str(restrictedend) + "\t" + querycoords + "/" + strand + "\t0.0\t" + strand + "\n"
                    aligncoveredregions = aligncoveredregions + coveredstring
                    coveredwindows = make_covered_window_string(ref, restrictedstart, restrictedend, query, querystart, queryend, strand, args)
        
                    aligncoveredwindows = aligncoveredwindows + coveredwindows
    
                    alignhetvariants = ""
                    novelhetvariantstring = ""
                    variantlist = alignparse.align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand)
                    for variant in variantlist:
                        chrom = variant.chrom
                        start = variant.start
                        end = variant.end
                        if start >= restrictedstart and end < restrictedend:
                            name = variant.name
                            alignment = alignname + "." + str(restrictedstart) + "_" + str(restrictedend)
                            variantstring = chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\t" + alignment + "\n"
                            alignhetvariants = alignhetvariants + variantstring
                            if name not in includedvariantdict.keys():
                                novelhetvariantstring = novelhetvariantstring + variantstring
                                includedvariantdict[name] = 1
    
                    hetvariants = hetvariants + novelhetvariantstring
    
                    variantintervals = pybedtools.BedTool(alignhetvariants, from_string=True)
                    coveredwindowintervals = pybedtools.BedTool(coveredwindows, from_string=True)
    
                    alignintervalswithcounts = bedtoolslib.intersectintervals(coveredwindowintervals, variantintervals, wa=True, counts=True)

                    numvars = len(variantintervals)
                    numwindows = len(coveredwindowintervals)
                    numcounts = len(alignintervalswithcounts)
                    logger.debug("Intersected " + str(numwindows) + " windows with " + str(numvars) + " variants and got " + str(numcounts) + " intervals with counts for align " + alignname)
                    windowintervalswithcounts = windowintervalswithcounts.cat(alignintervalswithcounts, postmerge=False)
                    logger.debug("windowintervalwithcounts has " + str(len(windowintervalswithcounts)) + " windows with counts for align " + alignname)

        logger.debug("Gathered " + str(len(windowintervalswithcounts)) + " windows with counts for align " + alignname)
        allwindowintervalswithcounts = allwindowintervalswithcounts.cat(windowintervalswithcounts, postmerge=False)

    return aligncoveredregions, allwindowintervalswithcounts, hetvariants

def make_covered_window_string(ref:str, refstart:int, refend:int, query:str, querystart:int, queryend:int, strand:str, args):
    windowsize = args.windowsize

    # if intervals is smaller than the windowsize, just return the interval as the window:
    if windowsize > refend - refstart:
        windowname = ref + "_" + str(refstart) + "_" + str(refend) + "_" + query + "_" + str(querystart) + "_" + str(queryend) + "_" + strand + "_0"
        bedstring = ref + "\t" + str(refstart) + "\t" + str(refend) + "\t" + windowname + "\n"
        return bedstring

    # otherwise return windows starting at "nice" positions with flanking partial windows:
    windowstart = windowsize*int(refstart/windowsize + 1.0)
    if windowstart - refstart > 0:
        bedstring = ref + "\t" + str(refstart) + "\t" + str(windowstart) + "\t" + ref + "_" + str(refstart) + "_" + str(refend) + "_" + query + "_" + str(querystart) + "_" + str(queryend) + "_" + strand + "_0" + "\n"
    else:
        bedstring = ""
    lastwindowend = windowsize*int(refend/windowsize)
    if refend < lastwindowend:
        lastwindowend = refend

    windownumber = 1
    while windowstart + windowsize <= lastwindowend:
        windowname = ref + "_" + str(refstart) + "_" + str(refend) + "_" + query + "_" + str(querystart) + "_" + str(queryend) + "_" + strand + "_" + str(windownumber)
        bedstring = bedstring + ref + "\t" + str(windowstart) + "\t" + str(windowstart + windowsize) + "\t" + windowname + "\n"
        windownumber = windownumber + 1
        windowstart = windowstart + windowsize

    if windowstart < refend:
        windowname = ref + "_" + str(refstart) + "_" + str(refend) + "_" + query + "_" + str(querystart) + "_" + str(queryend) + "_" + strand + "_" + str(windownumber)
        bedstring = bedstring + ref + "\t" + str(windowstart) + "\t" + str(refend) + "\t" + windowname + "\n"

    return bedstring

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    # read aligns from BAM files
    # coords dictionaries have keys that are the query coordinates and values that are the ref coordinates
    if args.bam1 is not None:
        alignobj1 = pysam.AlignmentFile(args.bam1, "rb")
        refbed1, querybed1, aligndict1 = read_aligns(alignobj1, args)
        querybed1.sort().saveas(args.prefix + ".hetquery1.bed")
    else:
        alignobj1 = None
    if args.bam2 is not None:
        alignobj2 = pysam.AlignmentFile(args.bam2, "rb")
        refbed2, querybed2, aligndict2 = read_aligns(alignobj2, args)
        querybed2.sort().saveas(args.prefix + ".hetquery2.bed")
    else:
        alignobj2 = None
    # refs
    if args.ref1 is not None:
        refobj1 = pysam.FastaFile(args.ref1)
    if args.ref2 is not None:
        refobj2 = pysam.FastaFile(args.ref2)

    # should we parse a gfa file? This is just an experimental option--not a true part of the program (random, I know)
    if args.gfa is not None:
        genomegraph = assemblygraph.read_graph_from_gfa(args.gfa, args)
        assemblygraph.find_hom_het_nodes(genomegraph, args)
        [compressed_ref1, chainfile1] = seqparse.compress_sequence(args.ref1, args)
        [compressed_ref2, chainfile2] = seqparse.compress_sequence(args.ref2, args)

        compare_compressed_fasta_to_node_seqs(compressed_ref1, args.gfa, args)
        compare_compressed_fasta_to_node_seqs(compressed_ref2, args.gfa, args)

    # look for corresponding alignments in both directions:
    if alignobj1 is not None:
        alignmapping = find_corresponding_alignment_pairs(alignobj1, alignobj2, refobj1, refobj2, aligndict1, aligndict2, querybed1, querybed2, args)
       
        coveredregions1, windowswithcounts1, hetvariants1 = find_hets_and_coveredregions(alignobj1, refobj1, refobj2, alignmapping, args)
        coveredregions2, windowswithcounts2, hetvariants2 = find_hets_and_coveredregions(alignobj2, refobj2, refobj1, alignmapping, args)

        logger.debug("find_hets_and_coveredregion returned " + str(len(windowswithcounts1)) + " and " + str(len(windowswithcounts2)) + " windows with counts")
        allvars = hetvariants1 + hetvariants2
        variantintervals = pybedtools.BedTool(allvars, from_string=True)
        variantintervalssorted = variantintervals.sort()
        #mergedvariants = bedtoolslib.mergeintervals(variantintervalssorted, collapsecolumns='4,5', collapseoutput='distinct', collapsedelim=',', distance=-1)
        variantintervalssorted.saveas(args.prefix + ".hetsites.bed")
        allcovered = coveredregions1 + coveredregions2
        coveredintervals = pybedtools.BedTool(allcovered, from_string=True)
        coveredintervals.sort().saveas(args.prefix + ".hapaligned.bed")
        windowwithcountintervals = windowswithcounts1.cat(windowswithcounts2, postmerge=False)
        windowwithcountintervals.sort().saveas(args.prefix + ".hapaligned.windows.withcounts.bed")
    
        if args.heatmap:
            heatmapwindows = heatmap.add_colors_to_intervals(windowwithcountintervals, args)
            if heatmapwindows is not None:
                heatmapwindows.sort().saveas(args.prefix + ".hapaligned.windows.withheatmap.bed")

if __name__ == "__main__":
    main()
