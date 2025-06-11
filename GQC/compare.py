import sys
import os
import re
import shutil
import pysam
from pysam import VariantFile
import argparse
import logging
import pybedtools
from pybedtools import BedTool
#import importlib.resources
from pathlib import Path
from GQC import output
from GQC import seqparse
from GQC import stats
from GQC import kmers
from GQC import align
from GQC import phasing
from GQC import bedtoolslib
from GQC import alignparse
from GQC import errors
from GQC import mummermethods
from GQC import structvar
from GQC import plots

logger = logging.getLogger(__name__)

def check_for_bedtools():
    if shutil.which("bedtools") is None:
        print("You don\'t seem to have bedtools in your path. Please install bedtools")
        logger.critical("You don\'t seem to have bedtools in your path. Please install bedtools")
        exit(1)
    return 0

def check_for_aligner(aligner:str):
    if aligner == "winnowmap2":
        aligner = "winnowmap"
    if shutil.which(aligner) is None:
        print("You don\'t seem to have the required aligner " + aligner + " installed in your path. Please install it")
        logger.critical("You don\'t seem to have the required aligner " + aligner + " installed in your path. Please install it")
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
        description="Compare two assemblies for the same diploid individual and report statistics"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('-Q', '--qname', type=str, required=False, default="AssemblyA", help='name of the assembly in which differences are reported (query)')
    parser.add_argument('-R', '--rname', type=str, required=False, default="AssemblyB", help='name of the assembly the first should be compared to (ref)')
    parser.add_argument('--q1fasta', type=str, required=True, help='path of a fasta file for haplotype 1 of the query assembly')
    parser.add_argument('--q2fasta', type=str, required=False, help='path of a fasta file for haplotype 2 of the query assembly')
    parser.add_argument('--r1fasta', type=str, required=True, help='path of a fasta file for haplotype 1 of the reference assembly')
    parser.add_argument('--r2fasta', type=str, required=False, help='path of a fasta file for haplotype 2 of the reference assembly')
    parser.add_argument('-p', '--prefix', type=str, required=True, help='prefix for output directory name, filenames use assembly names and this prefix (see -Q, -R)')
    parser.add_argument('-t', type=int, required=False, default=2, help='number of processors to use')
    parser.add_argument('-a', '--aligner', type=str, required=False, default='minimap2', help='aligner to use when comparing test assembly to reference assembly, can be minimap2 or winnowmap2 (default winnowmap2)')
    parser.add_argument('--alpha', type=float, required=False, default=0.05, help='emission probability for displaying haplotype markers for the wrong haplotype in HMM phase block algorithm')
    parser.add_argument('--beta', type=float, required=False, default=0.01, help='transition probability for changing haplotype state between adjacent markers (regardless of distance between them')
    parser.add_argument('-c', '--config', type=str, required=False, default="compareconfig.txt", help='path to a config file specifying locations of files used by GQC compare')
    parser.add_argument('--minns', type=int, required=False, default=10, help='minimum number of consecutive Ns required to break scaffolds into contigs')
    parser.add_argument('-m', '--minalignlength', type=int, required=False, default=500, help='minimum length of alignment required to be included in alignment statistics and error counts')
    parser.add_argument('--nosplit', action='store_true', required=False, help='use haplotype alignments without splitting on large indels. Default is to split alignments at locations with indels of at least --splitdistance')
    parser.add_argument('--splitdistance', type=int, required=False, default=10000, help='By default, split alignments when they contain indels of this size or greater')
    parser.add_argument('--maxclusterdistance', type=int, required=False, default=10000, help='maximum distance within a cluster of alignments')
    parser.add_argument('--vcf', action='store_true', required=False, default=False, help='write differences between assemblies in VCF (as well as BED) format')
    parser.add_argument('--haploid', action='store_true', required=False, help='run with just haploid assemblies q1 and r1')
    parser.add_argument('--hap2dip', action='store_true', required=False, help='compare a haploid assembly q1 to ref haplotypes r1 and r2')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    if not args.haploid and (not args.r1fasta or not args.r2fasta):
        logger.critical("Unless --haploid option is used, options r1fasta and r2fasta must be specified")
        exit(1)
    if not args.haploid and not args.hap2dip and not args.q2fasta:
        logger.critical("Unless --haploid or --hap2dip option is used, option q2fasta must be specified")
        exit(1)

    return args

def read_config_data(args)->dict:
    configfile = args.config
    configpath = Path(configfile)

    configvals = {}
    if configpath.exists():
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

    # log the command line:
    call_command = " ".join(sys.argv)
    logger.info(call_command)

    # check for necessary installed programs and write an output directory:
    check_for_bedtools()
    check_for_aligner(args.aligner)
    no_rscript = check_for_R()

    # dictionary of parameters from the configuration file (not used yet):
    compareparams = read_config_data(args)

    # pysam objects for the assembly haplotype fasta files:
    q1hapfasta = Path(args.q1fasta)
    r1hapfasta = Path(args.r1fasta)
    if (not q1hapfasta.is_file() or not r1hapfasta.is_file()):
        logger.critical("FASTA files specified with --q1fasta and --r1fasta must exist and be readable")
        exit(1)

    if not args.haploid:
        r2hapfasta = Path(args.r2fasta)
        if not r2hapfasta.is_file():
            logger.critical("FASTA file specified with --r2fasta must exist and be readable")
            exit(1)
   
    if not args.haploid and not args.hap2dip:
        q2hapfasta = Path(args.q2fasta)
        if not q2hapfasta.is_file():
            logger.critical("FASTA files specified with --q2fasta must exist and be readable")
            exit(1)
    
    hapq1obj = pysam.FastaFile(args.q1fasta)
    hapr1obj = pysam.FastaFile(args.r1fasta)
    hapdata = {'q1':{'pysamobj':hapq1obj, 'fasta':args.q1fasta, 'prefix':args.qname},
               'r1':{'pysamobj':hapr1obj, 'fasta':args.r1fasta, 'prefix':args.rname}}
    if not args.haploid:
        hapr2obj = pysam.FastaFile(args.r2fasta)
        hapdata['r1']['prefix'] = hapdata['r1']['prefix'] + "_hap1"
        hapdata['r2'] = {'pysamobj':hapr2obj, 'fasta':args.r2fasta, 'prefix':args.rname + "_hap2"}
        if not args.hap2dip:
            hapq2obj = pysam.FastaFile(args.q2fasta)
            hapdata['q1']['prefix'] = hapdata['q1']['prefix'] + "_hap1"
            hapdata['q2'] = {'pysamobj':hapq2obj, 'fasta':args.q2fasta, 'prefix':args.qname + "_hap2"}

    comparisondata = {'q1_to_r1':{'refobj':hapr1obj, 'queryobj':hapq1obj, 'refprefix':hapdata['r1']['prefix'], 'queryprefix':hapdata['q1']['prefix']}}
    if not args.haploid:
        comparisondata['q1_to_r2'] = {'refobj':hapr2obj, 'queryobj':hapq1obj, 'refprefix':hapdata['r2']['prefix'], 'queryprefix':hapdata['q1']['prefix']}
        if not args.hap2dip:
            comparisondata['q2_to_r1'] = {'refobj':hapr1obj, 'queryobj':hapq2obj, 'refprefix':hapdata['r1']['prefix'], 'queryprefix':hapdata['q2']['prefix']}
            comparisondata['q2_to_r2'] = {'refobj':hapr2obj, 'queryobj':hapq2obj, 'refprefix':hapdata['r2']['prefix'], 'queryprefix':hapdata['q2']['prefix']}

    comparisonoutputfiles = {}
    
    outputdir = output.create_output_directory(args.prefix)
    
    logger.info("Step 1 (of n): Scanning assembly sequences for N stretches, contig and scaffold lengths, etc.")
    bedregiondict = {}
    
    for haplotype in hapdata.keys():
        hapdict = hapdata[haplotype]
        seqparse.write_assembly_bedfiles(hapdict['pysamobj'], args, compareparams, hapdict['prefix'], bedregiondict)
        stats.write_assembly_haplotype_stats(hapdict['pysamobj'], bedregiondict[hapdict['prefix'] + "nonnregions"], bedregiondict[hapdict['prefix'] + "nregions"], args)
   

    q1_to_r1_phaseblockints = []
    q1phaseblockints = []
    q1_to_r2_phaseblockints = []
    q2_to_r1_phaseblockints = []
    q2_to_r2_phaseblockints = []

    if not args.haploid:
        logger.info("Step 2 (of n): Finding haplotype-specific kmers for each assembly")
    
        if not os.path.exists(outputdir + "/r1_not_r2.kmers.k40.ktab") or not os.path.exists(outputdir + "/r2_not_r1.kmers.k40.ktab"):
            for haplotype in hapdata.keys():
                hapdict = hapdata[haplotype]
                kmers.create_kmer_database(hapdict['fasta'], outputdir, hapdict['prefix'])
       
            kmers.find_anotb_kmers(hapdata['r1']['prefix'], hapdata['r2']['prefix'], outputdir, 'r1_not_r2') 
            kmers.find_anotb_kmers(hapdata['r2']['prefix'], hapdata['r1']['prefix'], outputdir, 'r2_not_r1') 

            # now can cleanup original databases:
            for haplotype in hapdata.keys():
                hapdict = hapdata[haplotype]
                kmers.remove_kmer_database(outputdir, hapdict['prefix'])
       
        logger.info("Step 3 (of n): Writing bed files of kmer locations of " + args.rname + " haplotype markers in the " + args.qname + " assembly")
        q1hapmerbed = kmers.map_kmer_markers_onto_fasta(hapdata['q1']['fasta'], [outputdir + "/" + 'r1_not_r2.kmers.k40', outputdir + "/" + 'r2_not_r1.kmers.k40'], outputdir)
        if not args.hap2dip:
            q2hapmerbed = kmers.map_kmer_markers_onto_fasta(hapdata['q2']['fasta'], [outputdir + "/" + 'r1_not_r2.kmers.k40', outputdir + "/" + 'r2_not_r1.kmers.k40'], outputdir)
    
        # use HMM algorithm to find matching phase blocks between assemblies
        logger.info("Step 4 (of n): Using HMM with emission probability " + str(args.alpha) + " and transition probability " + str(args.beta) + " to find matching phase blocks between assemblies")
        alpha = args.alpha
        transitionprob = args.beta
        q1phaseblockbed = outputdir + "/" + hapdata['q1']['prefix'] + ".hmmphasedscaffolds.bed"
        q1phaseblockmergedbed = q1phaseblockbed.replace('.bed', '.merged.bed')
        logger.info("Writing phase block bed file of haplotype kmers of " + args.rname + " present within the " + args.qname + " assembly")
        q1phaseblockints = phasing.find_hapmer_phase_blocks_with_hmm(q1hapmerbed, q1phaseblockbed, hapdata['q1']['pysamobj'], alpha, transitionprob, 0)
        if not os.path.exists(q1phaseblockmergedbed):
            q1phaseblockints.saveas(q1phaseblockmergedbed)
        q1_to_r1_phaseblockints = q1phaseblockints.filter(lambda x: x.name=="r1")
        q1_to_r2_phaseblockints = q1phaseblockints.filter(lambda x: x.name=="r2")

        if not args.hap2dip:
            q2phaseblockbed = outputdir + "/" + hapdata['q2']['prefix'] + ".hmmphasedscaffolds.bed"
            q2phaseblockmergedbed = q2phaseblockbed.replace('.bed', '.merged.bed')
            q2phaseblockints = phasing.find_hapmer_phase_blocks_with_hmm(q2hapmerbed, q2phaseblockbed, hapdata['q2']['pysamobj'], alpha, transitionprob, 0)
            if not os.path.exists(q2phaseblockmergedbed):
                q2phaseblockints.saveas(q2phaseblockmergedbed)
            q2_to_r1_phaseblockints = q2phaseblockints.filter(lambda x: x.name=="r1")
            q2_to_r2_phaseblockints = q2phaseblockints.filter(lambda x: x.name=="r2")

        ##stats.write_phase_block_stats(phaseblockints, outputfiles, benchmark_stats, args)
        #q1_to_r1_numintervals = len(q1_to_r1_phaseblockints)
        #print("End of assigned " + str(q1_to_r1_numintervals))

    else:
        logger.info("Skipping steps 2 through 4--not needed for comparing haploid assemblies")
   
    # align test assembly separately to each ref haplotype separately:
    logger.info("Step 5 (of n): Aligning test haplotypes separately to reference haplotypes and gathering trimmed alignments of phased assembly regions to their corresponding reference haplotype region")
    if not os.path.exists(hapdata['r1']['prefix'] + ".merge.sort.bam"):
        q1_to_r1_prefix = hapdata['q1']['prefix'] + "_to_" + hapdata['r1']['prefix'] + "." + args.aligner
        q1_to_r1_bamfile = align.align_haplotype_to_haplotype(hapdata['q1']['fasta'], hapdata['r1']['fasta'], q1_to_r1_prefix, compareparams, args)

        if not args.haploid:
            q1_to_r2_prefix = hapdata['q1']['prefix'] + "_to_" + hapdata['r2']['prefix'] + "." + args.aligner
            q1_to_r2_bamfile = align.align_haplotype_to_haplotype(hapdata['q1']['fasta'], hapdata['r2']['fasta'], q1_to_r2_prefix, compareparams, args)
            if not args.hap2dip:
                q2_to_r1_prefix = hapdata['q2']['prefix'] + "_to_" + hapdata['r1']['prefix'] + "." + args.aligner
                q2_to_r1_bamfile = align.align_haplotype_to_haplotype(hapdata['q2']['fasta'], hapdata['r1']['fasta'], q2_to_r1_prefix, compareparams, args)
                q2_to_r2_prefix = hapdata['q2']['prefix'] + "_to_" + hapdata['r2']['prefix'] + "." + args.aligner
                q2_to_r2_bamfile = align.align_haplotype_to_haplotype(hapdata['q2']['fasta'], hapdata['r2']['fasta'], q2_to_r2_prefix, compareparams, args)
    
            # read in alignments from BAM format, filtering out secondaries and finding the optimal alignment on the correct haplotype for each phase block in the assembly
            logger.info("Finding subaligns in query alignments for haplotype phase blocked regions of the assembly")
            q1_to_r1_trimmedbamfile = q1_to_r1_bamfile.replace(".bam", ".trimmed.bam")
            q1_to_r1_trimmedsortbamfile = q1_to_r1_bamfile.replace(".bam", ".trimmed.sort.bam")
            if not os.path.exists(q1_to_r1_trimmedsortbamfile):
                alignparse.trim_bamfile_to_intervals(q1_to_r1_bamfile, q1_to_r1_phaseblockints, q1_to_r1_trimmedbamfile, q1_to_r1_bamfile, args, sort=True, index=True)
            comparisondata['q1_to_r1']['bamfile'] = q1_to_r1_trimmedsortbamfile

            q1_to_r2_trimmedbamfile = q1_to_r2_bamfile.replace(".bam", ".trimmed.bam")
            q1_to_r2_trimmedsortbamfile = q1_to_r2_bamfile.replace(".bam", ".trimmed.sort.bam")
            if not os.path.exists(q1_to_r2_trimmedsortbamfile):
                alignparse.trim_bamfile_to_intervals(q1_to_r2_bamfile, q1_to_r2_phaseblockints, q1_to_r2_trimmedbamfile, q1_to_r2_bamfile, args, sort=True, index=True)
            comparisondata['q1_to_r2']['bamfile'] = q1_to_r2_trimmedsortbamfile

            if not args.hap2dip:
                q2_to_r1_trimmedbamfile = q2_to_r1_bamfile.replace(".bam", ".trimmed.bam")
                q2_to_r1_trimmedsortbamfile = q2_to_r1_bamfile.replace(".bam", ".trimmed.sort.bam")
                if not os.path.exists(q2_to_r1_trimmedsortbamfile):
                    alignparse.trim_bamfile_to_intervals(q2_to_r1_bamfile, q2_to_r1_phaseblockints, q2_to_r1_trimmedbamfile, q2_to_r1_bamfile, args, sort=True, index=True)
                comparisondata['q2_to_r1']['bamfile'] = q2_to_r1_trimmedsortbamfile
    
                q2_to_r2_trimmedbamfile = q2_to_r2_bamfile.replace(".bam", ".trimmed.bam")
                q2_to_r2_trimmedsortbamfile = q2_to_r2_bamfile.replace(".bam", ".trimmed.sort.bam")
                if not os.path.exists(q2_to_r2_trimmedsortbamfile):
                    alignparse.trim_bamfile_to_intervals(q2_to_r2_bamfile, q2_to_r2_phaseblockints, q2_to_r2_trimmedbamfile, q2_to_r2_bamfile, args, sort=True, index=True)
                comparisondata['q2_to_r2']['bamfile'] = q2_to_r2_trimmedsortbamfile
    
        else: 
            comparisondata['q1_to_r1']['bamfile'] = q1_to_r1_bamfile
 
    logger.info("Step 6 (of 11): Filtering alignments to include primary best increasing subset")

    for comparison in comparisondata.keys():
        trimmedphasedbam = comparisondata[comparison]['bamfile']
        refobj = comparisondata[comparison]['refobj']
        queryobj = comparisondata[comparison]['queryobj']
        if not args.nosplit:
            splitbam_name = trimmedphasedbam.replace(".bam", ".split.bam")
            splitsortbam_name = splitbam_name.replace(".bam", ".sort.bam")
            if not os.path.exists(splitsortbam_name):
                alignobj = pysam.AlignmentFile(trimmedphasedbam, "rb")
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
        logger.info("Step 5a (of 11): Assessing overall structural alignment of assembly")
        # (by default, rlis_aligndata are split alignments filtered for RLIS)
        outputfiles = {"alignplotdir":outputdir + "/alignmentplots", "alignplotprefix":outputdir + "/alignmentplots/" + comparison + ".clustered_aligns",
                "structvariantbed":outputdir + "/" + comparison + ".svs.bed"}
        if comparison not in comparisonoutputfiles.keys():
            comparisonoutputfiles[comparison] = {}
        comparisonoutputfiles[comparison]['alignplotprefix'] = outputdir + "/alignmentplots/" + comparison + ".clustered_aligns"
        comparisonoutputfiles[comparison]['structvariantbed'] = outputdir + "/" + comparison + ".svs.bed"
        bedregiondict["allexcludedregions"] = None
        benchmark_stats = {}
        if not os.path.exists(comparisonoutputfiles[comparison]['structvariantbed']):
            alignparse.assess_overall_structure(rlis_aligndata, refobj, queryobj, outputfiles, bedregiondict, benchmark_stats, args)
            structvar.write_structural_errors(rlis_aligndata, refobj, queryobj, outputfiles, benchmark_stats, args)
        else:
            logger.info("Not writing structural variant bed file " + comparisonoutputfiles[comparison]['structvariantbed'] + " because it already exists!")

        #stats.write_aligned_cluster_stats(outputfiles, benchmark_stats, args)

        # arguments only used in assembly benchmarking:
        pafaligns = None
        hetsites = None
        querycoveredbedfile = outputdir + "/" + comparison + ".querycovered.bed"
        refcoveredbedfile = outputdir + "/" + comparison + ".refcovered.bed"
        outputpat = None
        excludedregions = None
        hetarraybed = None

        [refcoveredbed, querycoveredbed, variants, hetsitealleles, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts] = alignparse.write_bedfiles(alignobj, pafaligns, refobj, queryobj, hetsites, querycoveredbedfile, outputpat, refcoveredbedfile, hetarraybed, excludedregions, args)

        # create merged unique outputfiles:
        [mergedrefcoveredbed, mergedrefcoveredbedfile] = bedtoolslib.mergebed(refcoveredbedfile)
        [mergedquerycoveredbed, mergedquerycoveredbedfile] = bedtoolslib.mergebed(querycoveredbedfile)

        comparisonoutputfiles[comparison]['refcoveredbed'] = refcoveredbed
        comparisonoutputfiles[comparison]['refcoveredbedfile'] = refcoveredbedfile
        comparisonoutputfiles[comparison]['querycoveredbed'] = querycoveredbed
        comparisonoutputfiles[comparison]['querycoveredbedfile'] = querycoveredbedfile
        comparisonoutputfiles[comparison]['mergedrefcoveredbed'] = mergedrefcoveredbed
        comparisonoutputfiles[comparison]['mergedquerycoveredbed'] = mergedquerycoveredbed

        logger.info("Step 7 (of 11): Writing primary alignment statistics about " + comparisondata[comparison]['queryprefix'] + " aligned to " + comparisondata[comparison]['refprefix'])
        #stats.write_merged_aligned_stats(refobj, queryobj, mergedtruthcoveredbed, mergedtestmatcoveredbed, mergedtestpatcoveredbed, outputfiles, benchmark_stats, args)

        if alignobj is not None:
            # classify variant errors as phasing or novel errors:
            logger.info("Step 8 (of 11): Writing discrepancy files")
            #stats.write_het_stats(outputfiles, benchmark_stats, args)
            outputfiles["testerrortypebed"] = outputdir + "/" + comparison + ".querydiscrepancies.bed"
            outputfiles["bencherrortypebed"] = outputdir + "/" + comparison + ".refdiscrepancies.bed"
            if not os.path.exists(outputfiles["testerrortypebed"]) or not os.path.exists(outputfiles["bencherrortypebed"]):
                errors.classify_errors(refobj, queryobj, variants, hetsites, outputfiles, compareparams, benchmark_stats, args)
            else:
                logger.info("Skipping writing discrepancy files--" + outputfiles["testerrortypebed"] + " and " + outputfiles["bencherrortypebed"] + " already exist!")
            comparisonoutputfiles[comparison]['queryerrorbed'] = outputdir + "/" + comparison + ".querydiscrepancies.bed"
            comparisonoutputfiles[comparison]['referrorbed'] = outputdir + "/" + comparison + ".refdiscrepancies.bed"
            if args.vcf:
                comparisonoutputfiles[comparison]['referrorvcf'] = outputdir + "/" + comparison + ".refdiscrepancies.vcf"
            #stats.write_qv_stats(benchmark_stats, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts, outputfiles, args)

    # Report general statistics across all haplotype comparisons:
    #stats.write_comparison_stats_file(hapdata, comparisondata, comparisonoutputfiles)

    # combine BED/VCF files:
    # Structural variants:
    combinedsvfile = outputdir + "/" + args.qname + "_vs_" + args.rname + ".svs.sort.bed"
    if not os.path.exists(combinedsvfile):
        combinedsvbedstring = ""
        for comparison in comparisondata.keys():
            with open(comparisonoutputfiles[comparison]['structvariantbed']) as svfh:
                combinedsvbedstring = combinedsvbedstring + svfh.read()
        combinedsvobj = pybedtools.bedtool.BedTool(combinedsvbedstring, from_string=True)
        combinedsvobj.sort().saveas(combinedsvfile)

    # Combined coverage of reference:
    combinedcovfile = outputdir + "/" + args.qname + "_vs_" + args.rname + ".refcovered.sort.bed"
    if not os.path.exists(combinedcovfile):
        combinedcovbedstring = ""
        for comparison in comparisondata.keys():
            with open(comparisonoutputfiles[comparison]['refcoveredbedfile']) as covfh:
                combinedcovbedstring = combinedcovbedstring + covfh.read()
        combinedcovobj = pybedtools.BedTool(combinedcovbedstring, from_string=True)
        combinedcovobj.sort().saveas(combinedcovfile)
    else:
        combinedcovobj = pybedtools.BedTool(combinedcovfile)

    # Combined coverage of query:
    combinedcovfile = outputdir + "/" + args.rname + "_vs_" + args.qname + ".querycovered.sort.bed"
    if not os.path.exists(combinedcovfile):
        combinedcovbedstring = ""
        for comparison in comparisondata.keys():
            with open(comparisonoutputfiles[comparison]['querycoveredbedfile']) as covfh:
                combinedcovbedstring = combinedcovbedstring + covfh.read()
        combinedcovobj = pybedtools.BedTool(combinedcovbedstring, from_string=True)
        combinedcovobj.sort().saveas(combinedcovfile)
    else:
        combinedcovobj = pybedtools.BedTool(combinedcovfile)

    # Uncovered regions in the reference:
    uncoveredfile = outputdir + "/" + args.qname + "_vs_" + args.rname + ".refuncovered.sort.bed"
    if not os.path.exists(uncoveredfile):
        hap1prefix = hapdata['r1']['prefix']
        hap2prefix = hapdata['r2']['prefix']
        hap1bedobj = bedregiondict[hap1prefix + "genomeregions"]
        hap2bedobj = bedregiondict[hap2prefix + "genomeregions"]
        genomebedtool = hap1bedobj.cat(hap2bedobj)

        uncoveredbedtool = bedtoolslib.subtractintervals(genomebedtool, combinedcovobj)
        uncoveredbedtool.sort().saveas(uncoveredfile)
    
    # Combined reference error bed files:
    combinedreferrorfile = outputdir + "/" + args.qname + "_vs_" + args.rname + ".refdiscrepancies.sort.bed"
    if not os.path.exists(combinedreferrorfile):
        combinedreferrorbedstring = ""
        for comparison in comparisondata.keys():
            with open(comparisonoutputfiles[comparison]['referrorbed']) as referrorfh:
                logger.info("Including ref discrepancy file" + comparisonoutputfiles[comparison]['referrorbed'])
                combinedreferrorbedstring = combinedreferrorbedstring + referrorfh.read()
        combinedreferrorobj = pybedtools.BedTool(combinedreferrorbedstring, from_string=True)
        logger.info("Saving to " + combinedreferrorfile)
        combinedreferrorobj.sort().saveas(combinedreferrorfile)
    else:
        combinedreferrorobj = pybedtools.BedTool(combinedreferrorfile)
    
    # Combined reference error VCF files:
    #if args.vcf:
        #combinedreferrorvcffile = outputdir + "/" + args.qname + "_vs_" + args.rname + ".refdiscrepancies.sort.vcf"
        #if not os.path.exists(combinedreferrorvcffile):
            #combinedheader = ""
            #headervcffile = comparisonoutputfiles['q1_to_r1']['referrorvcf']
            #headervcf_in = VariantFile(headervcffile)
            #combinedheader = headervcf_in.header
            #vcf_out = VariantFile(combinedreferrorvcffile, "w", header=combinedheader)
            #logger.info("Writing out to " + combinedreferrorvcffile)
            #for comparison in comparisondata.keys():
                #vcffile = comparisonoutputfiles[comparison]['referrorvcf']
                #logger.info("Reading in VCF " + vcffile)
                #vcf_in = VariantFile(vcffile)
                #for rec in vcf_in.fetch():
                    #vcf_out.write(rec)
            ##vcf_out.close()

            #with open(comparisonoutputfiles[comparison]['referrorvcf']) as referrorfh:
                #for line in referrorfh:
                    #if not startwwithpound:
                        #combinedreferrorvcfstring = combinedreferrorvcfstring + referrorfh.read()
        #combinedreferrorobj = pybedtools.BedTool(combinedreferrorbedstring, from_string=True)
        #logger.info("Saving to " + combinedreferrorfile)
        #combinedreferrorobj.sort().saveas(combinedreferrorfile)
    #else:
        #combinedreferrorobj = pybedtools.BedTool(combinedreferrorfile)
#
    ## plot alignment referrorerage across assembly and genome:
    if not no_rscript:
        #logger.info("Step 11 (of 11): Creating plots")
        #if not args.structureonly:
            #plots.plot_benchmark_align_coverage(args.assembly, args.benchmark, outputdir, benchparams)
            #plots.plot_testassembly_align_coverage(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
            #plots.plot_assembly_error_stats(args.assembly, args.benchmark, outputdir)
            #if alignobj is not None:
                #plots.plot_mononuc_accuracy(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
                #if len(alignedscorecounts) > 0:
                    #plots.plot_qv_score_concordance(args.assembly, args.benchmark, outputdir, benchparams["resourcedir"])
        for comparison in comparisondata.keys():
            refobj = comparisondata[comparison]['refobj']
            plots.plot_svcluster_align_plots(args.qname, args.rname, outputfiles["alignplotdir"], refobj, mode='compare', prefix=comparison)


if __name__ == "__main__":
    main()
