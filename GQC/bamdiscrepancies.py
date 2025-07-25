import sys
import os
import pysam
import argparse
import pybedtools
import logging
from pathlib import Path
from collections import namedtuple
from GQC import alignparse
from GQC import errors

logger = logging.getLogger(__name__)

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="For a bam file, report all within-alignment discrepancies (single base mismatches, insertions, and deletions"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('--bam', required=True, help='bam file of alignments of sequences mapped to the reference fasta specified with --ref')
    parser.add_argument('--ref', type=str, required=True, help='(indexed) reference fasta file for the bam file reference')
    parser.add_argument('--query', type=str, required=True, help='(indexed) query fasta file that was aligned to the bam file reference')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in het site gathering')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="bamdiscrepancies", help='prefix to use in output filenames')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

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

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output to ' + logfile + ' for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    # read aligns from the BAM file
    alignobj = pysam.AlignmentFile(args.bam, "rb")

    # reference
    refobj = pysam.FastaFile(args.ref)
    queryobj = pysam.FastaFile(args.query)

    alignmapping = {}
    variantnum = 1
    for align in alignobj.fetch():
        if align.is_secondary:
            continue
        else:
            print("Primary alignment")

        query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
        refcoords = ref + ":" + str(refstart) + "-" + str(refend)
        querycoords = query + ":" + str(querystart) + "-" + str(queryend)
        logger.info("Parsing alignment between " + querycoords + " and " + refcoords)
        allvars = alignparse.align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand)
        numvars = len(allvars)
        logger.info("Found " + str(numvars) + " variants in align with ref " + refcoords + " and query " + querycoords)
        for variant in allvars:
            #newchrom = variant.chrom
            #newstart = variant.start
            #newend = variant.end
            #newname = args.prefix + "." + str(variantnum)
            #newvartype = variant.vartype
            #newexcluded = variant.excluded
            #newqvscore = variant.qvscore
            #newvariant = varianttuple(chrom=newchrom, start=newstart, end=newend, name=newname, vartype=newvartype, excluded=newexcluded, qvscore=newqvscore)
            variantnum = variantnum + 1
            vcfrecord = errors.vcf_format(variant, refobj, queryobj)
            print(vcfrecord, end="")


if __name__ == "__main__":
    main()
