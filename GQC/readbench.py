import sys
import os
import re
import shutil
import pysam
import argparse
import pybedtools
import logging
import importlib.resources
from pathlib import Path
from GQC import errors
from GQC import output
from GQC import seqparse
from GQC import alignparse
from GQC import phasing
from GQC import stats
from GQC import plots

logger = logging.getLogger(__name__)

def check_for_R():
    if shutil.which("Rscript") is None:
        logger.warning("You don\'t seem to have Rscript in your path. Plots will not be generated")
        return 1
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print read quality statistics from a bam file of the reads aligned to a benchmark assembly"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-b', '--bam', required=True, default=None, help='bam file of alignments of the reads to the diploid benchmark')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='(indexed) fasta file for benchmark reference')
    parser.add_argument('--regions', '--regionbed', type=str, required=False, default='', help='bed file of benchmark regions to assess (default is to assess the whole diploid benchmark genome)')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in alignment statistics and error counts')
    parser.add_argument('--mincontiglength', type=int, required=False, default=500, help='minimum length for contig to be included in contig statistics')
    parser.add_argument('--maxstrreads', type=int, required=False, default=500, help='maximum number of reads to be evaluated for any one STR run in the STR run accuracy analysis. Use 0 to analyze all reads')
    parser.add_argument('--nostrs', action='store_true', required=False, help='skip analysis of short tandem repeat length accuracy')
    parser.add_argument('--nobaseerrors', action='store_true', required=False, help='skip analysis of base errors within reads')
    parser.add_argument('--downsample', type=restricted_float, required=False, default=None, help='fraction of read alignments to include in error reporting statistics calculations (must be a floating point number between 0 and 1)')
    parser.add_argument('--strs', action='store_true', required=False, help='analyse short tandem repeat accuracy')
    parser.add_argument('-e', '--errorfile', type=str, required=False, default='', help='preexisting file of read errors to report and plot stats for')
    parser.add_argument('-R', '--readsetname', type=str, required=False, default="test", help='name of the assembly being tested--should be query sequence in bam file')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging')

    return parser

# function to check that argument is within range from 0 to 1
def restricted_float(x):
    try:
        x = float(x)
    except ValueError:
        raise argparse.ArgumentTypeError("%r not a floating-point literal" % (x,))

    if x <= 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range (0.0, 1.0]"%(x,))
    return x

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    if not args.bam:
        logger.critical("Must specify a bam file with --bam")
        exit(1)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if not configpath.exists():
        logger.critical("A config file must exist in the default location benchconfig.txt or be specified with the --config option.")
        exit(1)
        #logger.info("Using resource locations from default config file")
        #template_res = importlib.resources.files("GQC").joinpath('benchconfig.txt')
        #with importlib.resources.as_file(template_res) as configfile:
            #with open(configfile, "r") as cr:
                #configline = cr.readline()
                #while configline:
                    #p = re.compile(r'^([^#\s]+):+\s+(\S+)$')
                    #match = p.match(configline)
                    #if match:
                        #key = match.group(1)
                        #value = match.group(2)
                        #configvals[key] = value
                    #configline = cr.readline()
    else:
        logger.info("Using resource locations from " + configfile)
        with open(configfile, "r") as cr:
            configline = cr.readline()
            while configline:
                p = re.compile(r'^([^#\s]+):+\s+(\S+)$')
                match = p.match(configline)
                if match:
                    key = match.group(1)
                    value = match.group(2)
                    configvals[key] = value
                configline = cr.readline()

    return configvals

def main() -> None:

    args = parse_arguments(sys.argv[1:])
    #check_for_bedtools()
    no_rscript = check_for_R()
    
    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    outputdir = output.create_output_directory(args.prefix)
    benchparams = read_config_data(args)

    hetsites = phasing.read_hetsites(benchparams["hetsitevariants"])

    alignobj = pysam.AlignmentFile(args.bam, "rb")
    refobj = pysam.FastaFile(args.reffasta)

    outputfiles = output.name_read_stats_files(args, outputdir)

    benchintervals = seqparse.write_included_bedfile(refobj, args, benchparams, outputfiles)

    if args.downsample is not None:
        logger.info("Downsampling to fraction " + str(args.downsample) + " of reads")


    # evaluate STR runs:
    #if not args.nomononucs:
        #logger.info("Assessing accuracy of mononucleotide runs")
        #logger.debug(benchparams["mononucruns"])
        #logger.debug(outputfiles["mononucstatsfile"])
        #if args.newflankrules:
            #mononucstats = errors.assess_mononuc_read_coverage_require_flank(alignobj, refobj, benchparams["mononucruns"], outputfiles, benchintervals, hetsites, args)
        #else:
            #mononucstats = errors.assess_mononuc_read_coverage(alignobj, benchparams["mononucruns"], outputfiles, benchintervals, hetsites, args)
        #stats.write_read_mononuc_stats(mononucstats, outputfiles, args)
        #plots.plot_read_mononuc_stats(args.readsetname, args.benchmark, outputdir)

    if not args.nostrs:
        logger.info("Assessing accuracy of short tandem repeats")
        print("Assessing accuracy of short tandem repeats")
    
        if "tetranucruns" in benchparams:
            logger.debug(benchparams["tetranucruns"])
            tetranucstats = errors.assess_str_read_coverage_require_flank('tetranuc', alignobj, refobj, benchparams["tetranucruns"], outputfiles, benchintervals, hetsites, args)
            stats.write_read_str_stats('tetranuc', tetranucstats, outputfiles, args)
        else:
            logger.info("No tetranucleotide bed file specified in configuration file(tetranucruns)")
        
        if "trinucruns" in benchparams:
            logger.debug(benchparams["trinucruns"])
            trinucstats = errors.assess_str_read_coverage_require_flank('trinuc', alignobj, refobj, benchparams["trinucruns"], outputfiles, benchintervals, hetsites, args)
            stats.write_read_str_stats('trinuc', trinucstats, outputfiles, args)
        else:
            logger.info("No trinucleotide bed file specified in configuration file(trinucruns)")
    
        if "dinucruns" in benchparams:
            logger.debug(benchparams["dinucruns"])
            dinucstats = errors.assess_str_read_coverage_require_flank('dinuc', alignobj, refobj, benchparams["dinucruns"], outputfiles, benchintervals, hetsites, args)
            stats.write_read_str_stats('dinuc', dinucstats, outputfiles, args)
        else:
            logger.info("No dinucleotide bed file specified in configuration file(dinucruns)")
    
        if "mononucruns" in benchparams:
            logger.debug(benchparams["mononucruns"])
            mononucstats = errors.assess_str_read_coverage_require_flank('mononuc', alignobj, refobj, benchparams["mononucruns"], outputfiles, benchintervals, hetsites, args)
            stats.write_read_str_stats('mononuc', mononucstats, outputfiles, args)
        else:
            logger.info("No mononucleotide bed file specified in configuration file(mononucruns)")

    # evaluate errors within read alignments:
    if not args.nobaseerrors:
        logger.info("Assessing errors within read alignments")
        logger.debug(outputfiles["readerrorfile"])
        errorstats = errors.assess_read_align_errors(alignobj, refobj, outputfiles["readerrorfile"], benchintervals, hetsites, args)
        stats.write_read_error_summary(errorstats, outputfiles)
        if len(errorstats["alignedqualscorecounts"]) > 0:
            plots.plot_read_error_stats(args.readsetname, args.benchmark, outputdir)
        #if len(errorstats["totalpositioncounts"]) > 0:
            #logger.info("Not plotting position counts yet!")
            #plots.plot_error_position_stats(args.readsetname, args.benchmark, outputdir)


if __name__ == "__main__":
    main()
