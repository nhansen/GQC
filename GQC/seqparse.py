import os
import re
import pysam
import pybedtools
import logging
from pathlib import Path
from GQC import bedtoolslib
from collections import namedtuple
import time

logger = logging.getLogger(__name__)

def write_genome_bedfiles(queryobj, refobj, args, benchparams, outputfiles, bedobjects):

    write_excluded_bedfile(refobj, args, benchparams, outputfiles, bedobjects)
    write_test_genome_bedfile(queryobj, args, outputfiles, bedobjects)
    find_all_ns(queryobj, args, outputfiles, bedobjects)

# write_excluded_bedfile function is written for assembly benchmarking only, not read benchmarking:
def write_excluded_bedfile(refobj, args, benchparams, outputfiles, bedobjects):

    allexcludedbedfiles = []
    if args.excludefile is not None:
        allexcludedbedfiles.append(args.excludefile)
    if args.includefile is not None:
        write_nonincluded_file(refobj, args.includefile, outputfiles["nonincludedbed"])
        allexcludedbedfiles.append(outputfiles["nonincludedfile"])
    if "excluderegions" in benchparams.keys():
        configexcluderegions = benchparams["excluderegions"]
        configexcludepath = Path(configexcluderegions)
        if configexcludepath.is_file():
            allexcludedbedfiles.append(configexcluderegions)
        else:
            logger.warning("Exclude file " + configexcluderegions + " does not exist so will not be used")
    if "nstretchregions" in benchparams.keys():
        nexcluderegions = benchparams["nstretchregions"]
        nexcludepath = Path(nexcluderegions)
        if nexcludepath.is_file():
            allexcludedbedfiles.append(nexcluderegions)
        else:
            logger.warning("Benchmark N file " + nexcluderegions + " does not exist so will not be used")
    if len(allexcludedbedfiles) == 0:
        bedobjects["allexcludedregions"] = pybedtools.BedTool("", from_string = True)
    elif len(allexcludedbedfiles) == 1:
        bedobjects["allexcludedregions"] = pybedtools.BedTool(allexcludedbedfiles[0])
    elif len(allexcludedbedfiles) > 1:
        logger.info("Merging " + str(allexcludedbedfiles) + " to create a combined excluded regions file")
        bedobjects["allexcludedregions"] = bedtoolslib.mergemultiplebedfiles(allexcludedbedfiles)
    bedobjects["allexcludedregions"].saveas(outputfiles["allexcludedbed"])

    return 0

# write_included_bedfile function is written for read benchmarking only, not assembly benchmarking:
def write_included_bedfile(refobj, args, benchparams, outputfiles):

    if args.regions != "":
        benchintervals = pybedtools.BedTool(args.regions)
        logger.info("Restricting evaluation to read sequence aligning to regions in " + args.regions)
    else:
        genomebedstring = ""
        for chrom in refobj.references:
            chromlength = refobj.get_reference_length(chrom)
            genomebedstring = genomebedstring + chrom + "\t0\t" + str(chromlength) + "\n"
        benchintervals = pybedtools.BedTool(genomebedstring, from_string = True)

    if "excluderegions" in benchparams.keys():
        configexcluderegions = benchparams["excluderegions"]
        configexcludepath = Path(configexcluderegions)
        if configexcludepath.is_file():
            excludeintervals = pybedtools.BedTool(configexcludepath)
            logger.info("Not including read bases aligned to regions excluded in " + configexcluderegions)
            benchintervals = benchintervals.subtract(excludeintervals)

    benchintervals.saveas(outputfiles["includedbedfile"])
    logger.info("Benchmark regions included in this analysis are written to file " + outputfiles["includedbedfile"])

    return benchintervals

def write_test_genome_bedfile(queryobj, args, outputfiles, bedobjects):

    genomebedstring = ""
    for scaffold in queryobj.references:
        scaffoldlength = queryobj.get_reference_length(scaffold)
        scaffstring = scaffold + "\t0\t" + str(scaffoldlength) + "\n"
        genomebedstring += scaffstring
    bedobjects["testgenomeregions"] = pybedtools.BedTool(genomebedstring, from_string = True)
    bedobjects["testgenomeregions"].saveas(outputfiles["testgenomebed"])

    return 0

def write_nonincluded_file(refobj, includedbed, nonincludedbed):

    genomebedstring = ""
    for refentry in refobj.references:
        reflength = refobj.get_reference_length(refentry)
        genomebedstring = genomebedstring + refentry + "\t0\t" + reflength + "\n"

    bedobjects["benchgenomeregions"] = pybedtools.BedTool(genomebedstring, from_string = True)
    bedobjects["includedregions"] = pybedtools.BedTool(includedbed)
    bedobjects["nonincludedregions"] = bedobjects["benchgenomeregions"].subtract(bedobjects["includedregions"])
    bedobjects["nonincludedregions"].saveas(nonincludedbed)

    return 0

def find_all_ns(queryobj, args, outputfiles, bedobjects)->list:

    # did the command-line arguments specify a pre-existing bedfile of N-stretch locations?
    user_n_file = args.n_bedfile

    if not user_n_file:
        p = re.compile("N+")

    gapbedstring = ""
    contigbedstring = ""
    n_interval_dict = {}
    # if an n interval file has been specified in the command line options, read it into gapbedstring, and subtract from genome to get contigbedstring:
    if user_n_file:
        nbedobj = BedTool(user_n_file)
        for n_interval in nbedobj:
            chrom = n_interval.chrom
            start = n_interval.start
            end = n_interval.end
            gapname = n_interval.name
            gapbedstring += chrom + "\t" + str(start) + "\t" + str(end) + "\t" + gapname + "\n"
            if chrom in n_interval_dict.keys():
                n_interval_dict[chrom].append(n_interval)
            else:
                n_interval_dict[chrom] = [n_interval]
        for ref in queryobj.references:
            if ref not in n_interval_dict.keys():
                n_interval_dict[ref] = []
            contignum = 1
            contigstart = 0
            for interval in n_interval_dict[ref]:
                interval_name = ref + "." + str(contignum)
                contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(interval.start) + "\t" + interval_name + "\n"
                contignum = contignum + 1
                contigstart = interval.end
            refend = queryobj.get_reference_length(ref)
            interval_name = ref + "." + str(contignum)
            contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(refend) + "\t" + interval_name + "\n"
    else:
        for ref in queryobj.references:
            contignum = 1
            contigstart = 0
            findstring = 'N' * args.minns
            chromseq = queryobj.fetch(ref).upper()
            refend = queryobj.get_reference_length(ref)
            start = chromseq.find(findstring)
            while start != -1:
                end = start + args.minns
                while end < refend and chromseq[end] == 'N':
                    end = end + 1
                if end - start >= args.minns:
                    gapname = "N." + ref + "." + str(contignum)
                    gapstring = ref + "\t" + str(start) + "\t" + str(end) + "\t" + gapname + "\n"
                    gapbedstring += gapstring
                    contigend = start
                    contigname = ref + "." + str(contignum)
                    contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(contigend) + "\t" + contigname + "\n"
                    contignum = contignum + 1
                    contigstart = end
                start = chromseq.find(findstring, end, refend)

            contigname = ref + "." + str(contignum)
            contigbedstring += ref + "\t" + str(contigstart) + "\t" + str(refend) + "\t" + contigname + "\n"

        bedobjects["testnonnregions"] = pybedtools.BedTool(contigbedstring, from_string = True)
        bedobjects["testnregions"] = pybedtools.BedTool(gapbedstring, from_string = True)

    if outputfiles["testnonnbed"]:
        bedobjects["testnonnregions"].saveas(outputfiles["testnonnbed"])
    if outputfiles["testnbed"]:
        bedobjects["testnregions"].saveas(outputfiles["testnbed"])

    return 0

def revcomp(seq:str) -> str:
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'M': 'N', 'K': 'N', 'R': 'N', 'W': 'N', 'Y': 'N'}
    bases = list(seq)
    bases = bases[::-1]
    bases = [complement[base.upper()] for base in bases]

    return ''.join(bases)

def compress_sequence(fastafile:str)->str:
    refobj = pysam.FastaFile(fastafile)

    if args.prefix is not None:
        chainfile = args.prefix + ".hpc.chain"
        compressedfasta = args.prefix + ".hpc.fasta"
    else:
        chainfile = re.sub("\\.fa.*$", ".hpc.chain", fastafile)
        chainfile = re.sub(".*/", "", chainfile)
        compressedfasta = re.sub("\\.fa.*$", ".hpc.fasta", fastafile)
        compressedfasta = re.sub(".*/", "", compressedfasta)

    if os.path.isfile(chainfile) and os.path.isfile(compressedfasta):
        print("Using pre-existing chain file " + chainfile + " and compressed fasta file " + compressedfasta)
    else:
        print("Writing chain file and homopolymer compressed assembly fasta")

        cfh = open(chainfile, "w")
        with open(compressedfasta, "w") as hfh:
            refseqs = refobj.references
            repeatnum = 0
            chainid = 0
            for entry in refseqs:
                print("Processing " + entry)
                fullseq = refobj.fetch(entry).upper()
                reflength = len(fullseq)
                matches = 0
                insertions = 0
                lastbase = ""
                inhpc = False
                compressedstring = ""
                chains = []
                for base in fullseq:
                    if base != lastbase:
                        if insertions > 0:
                            chains.append(str(matches) + "\t0\t" + str(insertions) + "\n")
                            insertions = 0
                            matches = 0
                        matches = matches + 1
                        compressedstring = compressedstring + base
                        lastbase = base
                    else: # hp!
                        insertions = insertions + 1
                if matches > 0 and insertions == 0:
                    chains.append(str(matches) + "\n\n")
                elif matches > 0 and insertions > 0:
                    chains.append(str(matches) + "\t0\t" + str(insertions) + "\n0\n\n")

                chainid = chainid + 1
                hpclength = len(compressedstring)
                cfh.write("chain 4000 " + entry + ".compressed " + str(hpclength) + " + 0 " + str(hpclength) + " " + entry + " " + str(reflength) + " + 0 " + str(reflength) + " " + str(chainid) + "\n")
                for chainstring in chains:
                    cfh.write(chainstring)
                hfh.write(">" + entry + ".compressed\n" + compressedstring + "\n")
        cfh.close()

    return [compressedfasta, chainfile]

