import sys
import os
import re
import shutil
import pysam
import argparse
import logging
from pathlib import Path
from GQC import output
from GQC import variants

logger = logging.getLogger(__name__)

def check_for_executable(programname:str):
    if shutil.which(programname) is None:
        print("You don\'t seem to have " + programname + " in your path. Please install " + programname)
        logger.critical("You don\'t seem to have " + programname + " in your path. Please install " + programname)
        exit(1)
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Construct a transformed genome fasta from a VCF file, masking regions not included in a high-confidence BED file"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-g', '--gvcf', required=False, default=None, help='gVCF file from which VCF variants and high confidence regions will be extracted')
    parser.add_argument('--mingq', type=int, required=False, default=10, help='if parsing a gVCF file, min GQ value to require to include variants and/or reference regions')
    parser.add_argument('--vcf', required=False, default=None, help='VCF file of variants to be applied to a reference when constructing the transformed genome')
    parser.add_argument('--hcregions', type=str, required=False, default=None, help='bed file of high confidence regions (in coordinates of the VCF reference)')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='(indexed) fasta file for the VCF file reference genome')
    parser.add_argument('--hap1chroms', type=str, required=False, default='chrY', help='comma-delimited list of reference entries that are only on haplotype 1')
    parser.add_argument('--hap2chroms', type=str, required=False, default='chrX', help='comma-delimited list of reference entries that are only on haplotype 2')
    parser.add_argument('-s', '--samplename', type=str, required=True, help='name of sample in VCF file')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name')
    parser.add_argument('-t', type=int, required=False, default=2, help='number of processors to use when running multi-threaded programs')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    # check for necessary installed programs and write an output directory:
    check_for_executable("bedtools")
    check_for_executable("bcftools")
    check_for_executable("liftOver")

    # did the user pass a single gVCF file or a VCF/high conf bed combo?
    if args.gvcf:
        gvcffile = args.gvcf
        logger.info("Creating VCF and high-confidence BED file from gVCF file " + gvcffile)
        # here we'll call a new routine to do this
        [vcffile, hcbedfile] = variants.create_vcf_and_hcbed_from_gvcf(gvcffile, args)
    elif args.vcf and args.hcregions:
        vcffile = Path(args.vcf)
        hcbedfile = Path(args.hcregions)
        if not vcffile.is_file():
            logger.critical("VCF file " + args.vcffile + " must exist and be readable")
            print("VCF file " + args.vcffile + " must exist and be readable", file=sys.stderr)
            exit(1)
    else:
        logger.critical("Either a gvcf file (--gvcf option) or both a VCF file and a high-confidence BED file (options --vcf, --hcregions) must be specified")
        print("VCF file " + args.vcffile + " must exist and be readable", file=sys.stderr)
        exit(1)


    # pysam objects for the benchmark and test assembly fasta files:
    ref = Path(args.reffasta)
    if not ref.is_file():
        logger.critical("Ref fasta file " + args.reffasta + " must exist and be readable")
        print("Ref fasta file " + args.reffasta + " must exist and be readable", file=sys.stderr)
        exit(1)
    outputdir = output.create_output_directory(args.prefix)

    # dictionary of this run's output file names:
    outputfiles = output.name_constructgenome_output_files(args, outputdir)

    # apply VCF variants to the reference, then filter non-high-confidence regions with Ns:
    variants.create_constructed_genome(args.reffasta, args.vcf, args.hcregions, outputdir, outputfiles, args)

if __name__ == "__main__":
    main()
