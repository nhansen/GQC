import sys
import os
import re
import shutil
import pysam
import argparse
import logging
from pybedtools import BedTool
import importlib.resources
from pathlib import Path
from collections import namedtuple
from GQC import bedtoolslib
from GQC import errors
from GQC import output
from GQC import seqparse
from GQC import align
from GQC import alignparse
from GQC import structvar
from GQC import phasing
from GQC import stats
from GQC import mummermethods
from GQC import plots

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

logger = logging.getLogger(__name__)

def check_for_bedtools():
    if shutil.which("bedtools") is None:
        print("You don\'t seem to have bedtools in your path. Please install bedtools")
        logger.critical("You don\'t seem to have bedtools in your path. Please install bedtools")
        exit(1)
    return 0

def check_for_R():
    if shutil.which("Rscript") is None:
        print("You don\'t seem to have Rscript in your path. Plots will not be generated")
        logger.warning("You don\'t seem to have Rscript in your path. Plots will not be generated")
        return 1
    return 0

def check_for_fastk():
    if shutil.which("FastK") is None or shutil.which("KmerMap") is None:
        print("You don\'t seem to have FastK with the KmerMap command in your path. These are necessary so scaffold phase blocks can be assessed")
        logger.critical("You don\'t seem to have FastK with the KmerMap command in your path. These are necessary so scaffold phase blocks can be assessed")
        exit(1)
    return 0

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Print assembly statistics from bam files of the assembly aligned to a benchmark assembly"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-b', '--bam', required=False, default=None, help='bam file of alignments of the test (haploid) assembly to the diploid benchmark')
    parser.add_argument('--paf', required=False, default=None, help='paf-formatted file of alignments of the test (haploid) assembly to the diploid benchmark')
    parser.add_argument('-r', '--reffasta', type=str, required=True, help='(indexed) fasta file for benchmark reference')
    parser.add_argument('-q', '--queryfasta', type=str, required=True, help='(indexed) fasta file for haploid or diploid test assembly')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name, filenames have assembly name in prefix (see -A)')
    parser.add_argument('-t', type=int, required=False, default=2, help='number of processors to use')
    parser.add_argument('-a', '--aligner', type=str, required=False, default='winnowmap2', help='aligner to use when comparing assembly to benchmark, can be minimap2 or winnowmap2 (default winnowmap2)')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=5000, help='minimum length of alignment required to be included in alignment statistics and error counts')
    parser.add_argument('--mincontiglength', type=int, required=False, default=500, help='minimum length for contig to be included in contig statistics')
    parser.add_argument('--minns', type=int, required=False, default=10, help='minimum number of consecutive Ns required to break scaffolds into contigs')
    parser.add_argument('--maxclusterdistance', type=int, required=False, default=10000, help='maximum distance within a cluster of alignments')
    parser.add_argument('--merquryblocks', action='store_true', required=False, help='flag option for calculating phase blocks of the test assembly scaffolds using the algorithm used in Merqury rather than an HMM')
    parser.add_argument('--shortnum', type=int, required=False, default=100, help='maximum number of short-range phase switches to allow in phase blocks (if --merquryblocks option is specified)')
    parser.add_argument('--shortlimit', type=int, required=False, default=20000, help='maximum lengths of short-range switches allowable in phase blocks (if --merquryblocks option is specified)')
    parser.add_argument('--includefile', type=str, required=False, default=None, help='bed file of benchmark locations to include in the evaluation (stretches of 10 or more Ns and regions excluded in the exclude file will still not be considered)')
    parser.add_argument('--excludefile', type=str, required=False, default=None, help='bed file of benchmark locations to exclude from consideration (in addition to stretches of 10 or more Ns and regions in the exclude file specified in the config file)')
    parser.add_argument('--vcf', action='store_true', required=False, default=False, help='write differences from benchmark in VCF format')
    parser.add_argument('-n', '--n_bedfile', type=str, required=False, default=None, help='pre-existing bedfile of locations of N-stretches splitting scaffolds into contigs')
    parser.add_argument('--variantfile', type=str, required=False, default=None, help='pre-existing file of variant locations in assembly compared to benchmark')
    parser.add_argument('--structureonly', action='store_true', required=False, help='analyse only the long-range structure of the assembly')
    parser.add_argument('-A', '--assembly', type=str, required=False, default="test", help='name of the assembly being tested--should correspond to query sequence in bam file and will be used in output file names')
    parser.add_argument('-B', '--benchmark', type=str, required=False, default="truth", help='name of the assembly being used as a benchmark--should be the reference sequence in the bam file')
    parser.add_argument('-c', '--config', type=str, required=False, default="benchconfig.txt", help='path to a config file specifying locations of benchmark data files')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')
    parser.add_argument('--nosplit', action='store_true', required=False, help='use haplotype alignments without splitting on large indels. Default is to split alignments at locations with indels of at least --splitdistance')
    parser.add_argument('--splitdistance', type=int, required=False, default=10000, help='By default, split alignments when they contain indels of this size or greater')
    parser.add_argument('--alpha', type=float, required=False, default=0.05, help='emission probability for displaying opposite haplotype markers in HMM phase block algorithm')
    parser.add_argument('--beta', type=float, required=False, default=0.01, help='transition probability for changing haplotype state between adjacent markers (regardless of distance between them, unless --distancemultiplier value is set) in HMM phase block algorithm')
    #parser.add_argument('--distancemultiplier', type=int, required=False, default=0, help='set the transition probability between two different states to the value of beta times the number of bases between the two markers divided by this value')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if not configpath.exists():
        logger.critical("A config file must exist in the default location benchconfig.txt or be specified with the --config option.")
        exit(1)
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

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output to ' + logfile + ' for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    # check for necessary installed programs and write an output directory:
    check_for_bedtools()
    no_rscript = check_for_R()

    # dictionary of parameters from the benchmark configuration file:
    benchparams = read_config_data(args)

    # pysam objects for the benchmark and test assembly fasta files:
    ref = Path(args.reffasta)
    query = Path(args.queryfasta)
    if not ref.is_file() or not query.is_file():
        logger.critical("Ref fasta file " + args.reffasta + " and query fasta file " + args.queryfasta + " must exist and be readable")
        exit(1)
    refobj = pysam.FastaFile(args.reffasta)
    queryobj = pysam.FastaFile(args.queryfasta)

    outputdir = output.create_output_directory(args.prefix)

    # dictionary of this run's output file names:
    outputfiles = output.name_output_files(args, outputdir)

    logger.info("Step 1 (of 11): Writing bed files for excluded regions, test assembly scaffold spans, lengths, N stretches, and contigs (ignoring stretches of less than " + str(args.minns) + " Ns)")
    bedregiondict = {}
    # merged excluded regions are saved as BedTool object "allexcludedregions" in bedregiondict here:
    seqparse.write_genome_bedfiles(queryobj, refobj, args, benchparams, outputfiles, bedregiondict)
   

    # find general stats about contig/scaffold lengths, N/L50's, etc.:
    logger.info("Step 2 (of 11): Writing general statistics about " + args.assembly + " assembly")
    benchmark_stats = stats.write_general_assembly_stats(refobj, queryobj, bedregiondict["testnonnregions"], bedregiondict["testnregions"], outputfiles, args)

    logger.info("Step 3 (of 11): Phasing assembly regions using benchmark's haplotype-specific kmers")
    if args.merquryblocks:
        markerbed = outputfiles["phasemarkerbed"]
    else:
        markerbed = outputfiles["mergedphasemarkerbed"]
    if not os.path.exists(markerbed):
        print("File " + markerbed + " does not exist--calling FASTK to generate one")
        logger.info("File " + markerbed + " does not exist--calling FASTK to generate one")
        check_for_fastk()
        phasing.map_benchmark_hapmers_onto_assembly(args.queryfasta, benchparams["matmarkerdb"], benchparams["patmarkerdb"], outputdir, outputfiles)

    if args.merquryblocks: # use Merqury phase block algorithm
        logger.info("Using Merqury algorithm with shortnum " + str(args.shortnum) + " and shortlimit " + str(args.shortlimit) + " to find phase blocks")
        shortnum = args.shortnum
        shortlimit = args.shortlimit
        phaseblockints = phasing.find_phase_blocks_from_marker_bed(markerbed, queryobj.references, shortnum, shortlimit)
    else: # use HMM algorithm to find phase blocks
        logger.info("Using HMM with emission probability " + str(args.alpha) + " and transition probability " + str(args.beta) + " to find phase blocks")
        alpha = args.alpha
        transitionprob = args.beta
        phaseblockints = phasing.find_hapmer_phase_blocks_with_hmm(markerbed, outputfiles["hmmphaseblockbed"], queryobj.references, alpha, transitionprob, 0)
    phaseblockints.saveas(outputfiles["phaseblockbed"])
    matphaseblockints = phaseblockints.filter(lambda x: x.name=="mat")
    patphaseblockints = phaseblockints.filter(lambda x: x.name=="pat")
    #stats.write_phase_block_stats(phaseblockints, outputfiles, benchmark_stats, args)

    # align test assembly separately to maternal and paternal haplotypes:
    logger.info("Step 4 (of 11): Aligning assembly separately to maternal and paternal haplotypes of the benchmark and gathering trimmed alignments of phased assembly regions to their corresponding haplotypes")
    if not os.path.exists(outputfiles["trimmedphasedalignprefix"] + ".merge.sort.bam"):
        [matbenchbamfile, patbenchbamfile] = align.align_assembly_to_benchmark_haplotypes(args.queryfasta, outputfiles, benchparams, args)
    
        # read in alignments from BAM format, filtering out secondaries and finding the optimal alignment on the correct haplotype for each phase block in the assembly

        benchdiploidheaderfile = benchparams["benchdiploidheader"]
    
        [mataligns, matalignedintervals] = alignparse.index_aligns_by_boundaries(matbenchbamfile, args)
        [pataligns, patalignedintervals] = alignparse.index_aligns_by_boundaries(patbenchbamfile, args)
    
        matalignedintervals.saveas(outputfiles["matalignedregions"])
        patalignedintervals.saveas(outputfiles["patalignedregions"])
    
        print("Finding subaligns in maternal alignments for maternal phase blocked regions of the assembly")
        matblocksubaligns = alignparse.find_phaseblock_subaligns(matphaseblockints, matalignedintervals, mataligns)
        mattrimmedbamfile = outputfiles["trimmedphasedalignprefix"] + ".mat.bam"
        alignparse.write_aligns_to_bamfile(mattrimmedbamfile, matblocksubaligns, headerbam=matbenchbamfile, sort=False)
        print("Finding subaligns in paternal alignments for paternal phase blocked regions of the assembly")
        patblocksubaligns = alignparse.find_phaseblock_subaligns(patphaseblockints, patalignedintervals, pataligns)
        pattrimmedbamfile = outputfiles["trimmedphasedalignprefix"] + ".pat.bam"
        alignparse.write_aligns_to_bamfile(pattrimmedbamfile, patblocksubaligns, headerbam=patbenchbamfile, sort=False)

        # now merge the maternal and paternal trimmed files to a single file with a diploid header, sort, and index:
        alignparse.merge_trimmed_bamfiles(mattrimmedbamfile, pattrimmedbamfile, benchdiploidheaderfile, outputfiles)
    else: 
        logger.info("Skipping step 4 (of 11): Trimmed phased alignments already exist in " + outputfiles["trimmedphasedalignprefix"] + ".merge.sort.bam")

    trimmedphasedbam = outputfiles["trimmedphasedalignprefix"] + ".merge.sort.bam"

    logger.info("Step 5 (of 11): Filtering alignment to include primary best increasing subset")
    if not args.nosplit:
        alignobj = pysam.AlignmentFile(trimmedphasedbam, "rb")
        splitbam_name = trimmedphasedbam.replace(".bam", ".split.bam")
        splitsortbam_name = splitbam_name.replace(".bam", ".sort.bam")
        alignparse.split_aligns_and_sort(splitbam_name, alignobj, minindelsize=args.splitdistance)
        pysam.sort("-o", splitsortbam_name, splitbam_name)
        pysam.index(splitsortbam_name)
        alignobj = pysam.AlignmentFile(splitsortbam_name, "rb")
        aligndata = alignparse.read_bam_aligns(alignobj, args.minalignlength)
        rlis_aligndata = mummermethods.filter_aligns(aligndata, "target")
    else:
        alignobj = pysam.AlignmentFile(bamfile, "rb")
        aligndata = alignparse.read_bam_aligns(alignobj, args.minalignlength)
        rlis_aligndata = mummermethods.filter_aligns(aligndata, "target")

    ## find clusters of consistent, covering alignments and calculate continuity statistics:
    logger.info("Step 5 (of 11): Assessing overall structural alignment of assembly")
    # (by default, rlis_aligndata are split alignments filtered for RLIS)
    alignparse.assess_overall_structure(rlis_aligndata, refobj, queryobj, outputfiles, bedregiondict, benchmark_stats, args)
    structvar.write_structural_errors(rlis_aligndata, refobj, queryobj, outputfiles, benchmark_stats, args)
    stats.write_aligned_cluster_stats(outputfiles, benchmark_stats, args)

    if not args.structureonly:
       logger.info("Step 6 (of 11): Writing bed files of regions covered by alignments of " + args.assembly + " to " + args.benchmark)
       ## read in locations of het variants in the benchmark:
       hetsites = phasing.read_hetsites(benchparams["hetsitevariants"])
       hetarrays = phasing.sort_chrom_hetsite_arrays(hetsites)

       pafaligns = None
       [refcoveredbed, querycoveredbed, variants, hetsitealleles, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts] = alignparse.write_bedfiles(alignobj, pafaligns, refobj, queryobj, hetarrays, outputfiles["testmatcovered"], outputfiles["testpatcovered"], outputfiles["truthcovered"], outputfiles["coveredhetsitealleles"], bedregiondict["allexcludedregions"], args)

       ## create merged unique outputfiles:
       [mergedtruthcoveredbed, outputfiles["mergedtruthcovered"]] = bedtoolslib.mergebed(outputfiles["truthcovered"])
       [mergedtestmatcoveredbed, outputfiles["mergedtestmatcovered"]] = bedtoolslib.mergebed(outputfiles["testmatcovered"])
       [mergedtestpatcoveredbed, outputfiles["mergedtestpatcovered"]] = bedtoolslib.mergebed(outputfiles["testpatcovered"])

       logger.info("Step 7 (of 11): Writing primary alignment statistics about " + args.assembly + " assembly")
       stats.write_merged_aligned_stats(refobj, queryobj, mergedtruthcoveredbed, mergedtestmatcoveredbed, mergedtestpatcoveredbed, outputfiles, benchmark_stats, args)

       if alignobj is not None:
           ## classify variant errors as phasing or novel errors:
           logger.info("Step 8 (of 11): Writing phase switch statistics")
           stats.write_het_stats(outputfiles, benchmark_stats, args)
           logger.info("Step 9 (of 11): Determining whether errors are switched haplotype or novel")
           errors.classify_errors(refobj, queryobj, variants, hetsites, outputfiles, benchparams, benchmark_stats, args)
           stats.write_qv_stats(benchmark_stats, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts, outputfiles, args)
#
           ## evaluate mononucleotide runs:
           logger.info("Step 10 (of 11): Assessing accuracy of mononucleotide runs")
           bedtoolslib.intersectbed(benchparams["mononucruns"], outputfiles["mergedtruthcovered"], outputfile=outputfiles["coveredmononucsfile"], writefirst=True)
           mononucswithvariantsbedfile = bedtoolslib.intersectbed(outputfiles["coveredmononucsfile"], outputfiles["bencherrortypebed"], outputfiles["mononucswithvariantsfile"], outerjoin=True, writeboth=True)
           mononucstats = errors.gather_mononuc_stats(outputfiles["mononucswithvariantsfile"], outputfiles["mononucstatsfile"])
           stats.write_mononuc_stats(mononucstats, outputfiles, benchmark_stats, args)

    # plot alignment coverage across assembly and genome:
    if not no_rscript:
        logger.info("Step 11 (of 11): Creating plots")
        if not args.structureonly:
            plots.plot_benchmark_align_coverage(args.assembly, args.benchmark, outputdir, benchparams)
            plots.plot_testassembly_align_coverage(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
            plots.plot_assembly_error_stats(args.assembly, args.benchmark, outputdir)
            if alignobj is not None:
                plots.plot_mononuc_accuracy(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
                if len(alignedscorecounts) > 0:
                    plots.plot_qv_score_concordance(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
        plots.plot_svcluster_align_plots(args.assembly, args.benchmark, outputfiles["alignplotdir"], benchparams["resourcedir"], refobj)


if __name__ == "__main__":
    main()
