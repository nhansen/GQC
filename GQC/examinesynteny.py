import sys
import os
import re
import shutil
import pysam
import argparse
import pybedtools
import importlib.resources
from pathlib import Path
from GQC import alignparse

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Examine bam files aligning one diploid haplotype to another and report on synteny/heterozygosity"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('--bam1', required=True, help='bam file of alignments in one direction from one haplotype mapped to the other')
    parser.add_argument('--bam2', required=True, help='bam file of alignments in the other direction')
    parser.add_argument('--ref1', type=str, required=True, help='(indexed) reference fasta file for the first bam file')
    parser.add_argument('--ref2', type=str, required=True, help='(indexed) reference fasta file for the first bam file')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in het site gathering')
    parser.add_argument('--hetwindowsize', type=int, required=False, default=20000, help='window size for reporting heterozygosity levels')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="hetsites", help='prefix to use in output filenames')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def read_aligns(bamobj, args):
    align_coords = {}

    for align in bamobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= args.minalignlength:
            query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
            querycoords = query + ":" + str(querystart) + "-" + str(queryend)
            refcoords = ref + ":" + str(refstart) + "-" + str(refend)
            querychrom = query.split("_")[0]
            refchrom = ref.split("_")[0]

            if querychrom == refchrom:
                align_coords[querycoords] = refcoords

    return align_coords

def find_hets_and_homregions(bamobj, refobj, queryobj, alignmap:dict, args):

    # strings which will be used to create bed files with pybedtools:
    aligncoveredregions = ""
    aligncoveredwindows = ""
    hetvariants = ""
    homregions = ""
    for align in bamobj.fetch():
        query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
        querycoords = query + ":" + str(querystart) + "-" + str(queryend)
        refcoords = ref + ":" + str(refstart) + "-" + str(refend)

        if querycoords in alignmap.keys() and alignmap[querycoords] == refcoords:
            # add query\tstart\tstop\tname\tscore\tstrand to covered string:
            coveredstring = ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querycoords + "/" + strand + "\t0.0\t" + strand + "\n"
            aligncoveredregions = aligncoveredregions + coveredstring
            coveredwindows = make_covered_window_string(ref, refstart, refend, query, querystart, queryend, strand, args)
            aligncoveredwindows = aligncoveredwindows + coveredwindows
            variantlist = alignparse.align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand)
            lastquerypos = querystart
            for variant in variantlist:
                chrom = variant.chrom
                start = variant.start
                end = variant.end
                name = variant.name
                variantstring = chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\n"
                hetvariants = hetvariants + variantstring

    return aligncoveredregions, aligncoveredwindows, hetvariants, homregions

def make_covered_window_string(ref:str, refstart:int, refend:int, query:str, querystart:int, queryend:int, strand:str, args):
    bedstring = ""
    windowsize = args.hetwindowsize

    windowstart = windowsize*int(refstart/windowsize + 1.0)
    lastwindowend = windowsize*int(refend/windowsize)
    windownumber = 1
    while windowstart < lastwindowend:
        windowname = ref + "_" + str(refstart) + "_" + str(refend) + "_" + query + "_" + str(querystart) + "_" + str(queryend) + "_" + strand + "_" + str(windownumber)
        bedstring = bedstring + ref + "\t" + str(windowstart) + "\t" + str(windowstart + windowsize) + "\t" + windowname + "\n"
        windownumber = windownumber + 1
        windowstart = windowstart + windowsize
    return bedstring

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    # read aligns from BAM files
    # coords dictionaries have keys that are the query coordinates and values that are the ref coordinates
    alignobj1 = pysam.AlignmentFile(args.bam1, "rb")
    coords1 = read_aligns(alignobj1, args)
    alignobj2 = pysam.AlignmentFile(args.bam2, "rb")
    coords2 = read_aligns(alignobj2, args)
    # refs
    refobj1 = pysam.FastaFile(args.ref1)
    refobj2 = pysam.FastaFile(args.ref2)

    alignmapping = {}
    for align1 in alignobj1.fetch():
        ref1, refstart1, refend1, query1, querystart1, queryend1, strand1 = alignparse.retrieve_align_data(align1)
        refcoords1 = ref1 + ":" + str(refstart1) + "-" + str(refend1)
        querycoords1 = query1 + ":" + str(querystart1) + "-" + str(queryend1)

        if querycoords1 in coords2.keys():
            refcoords2 = coords2[querycoords1]
            if refcoords1 == refcoords2 and refend1 - refstart1 >= 100000:
                alignmapping[refcoords1] = querycoords1
                alignmapping[querycoords1] = refcoords1

    coveredregions1, coveredwindows1, hetvariants1, homregions1 = find_hets_and_homregions(alignobj1, refobj1, refobj2, alignmapping, args)
    coveredregions2, coveredwindows2, hetvariants2, homregions2 = find_hets_and_homregions(alignobj2, refobj2, refobj1, alignmapping, args)

    allvars = hetvariants1 + hetvariants2
    variantintervals = pybedtools.BedTool(allvars, from_string=True)
    variantintervals.sort().saveas(args.prefix + ".hetsites.bed")
    allcovered = coveredregions1 + coveredregions2
    coveredintervals = pybedtools.BedTool(allcovered, from_string=True)
    coveredintervals.sort().saveas(args.prefix + ".1to1alignable.bed")
    allcoveredwindows = coveredwindows1 + coveredwindows2
    coveredwindowintervals = pybedtools.BedTool(allcoveredwindows, from_string=True)
    coveredwindowintervals.sort().saveas(args.prefix + ".1to1alignable.windows.bed")

if __name__ == "__main__":
    main()
