import os
import sys
import subprocess
import re
import pysam
import pybedtools
import logging
import math
from pathlib import Path
from collections import namedtuple
from pybedtools import BedTool
from GQC import seqparse
from GQC import kmers
from GQC import bedtoolslib

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

logger = logging.getLogger(__name__)

def read_hetsites(hetsitefile)->dict:
    hets = Path(hetsitefile)
    hetsites = {}
    if not hets.is_file():
        logger.error("Het site file " + hetsitefile + " does not exist so phasing analysis will be incomplete")
    else:
        with open(hetsitefile, "r") as hfh:
            hetsiteline = hfh.readline()
            while hetsiteline:
                hetsiteline = hetsiteline.rstrip()
                [chrom, start, end, name] = hetsiteline.split("\t")
                namefields = name.split("_")
                refallele = namefields[-3]
                altallele = namefields[-2]
                hetsitename = chrom + "_" + str(int(start) + 1) + "_" + refallele + "_" + altallele 
                if refallele != "*" and altallele != "*" and len(refallele)==1 and len(altallele)==1:
                    vartype = 'SNV'
                else:
                    vartype = 'INDEL'
                hetsites[hetsitename] = varianttuple(chrom=chrom, start=int(start), end=int(end), name=name, vartype=vartype, excluded=False, qvscore=None)
                hetsiteline = hfh.readline()

    return hetsites

def sort_chrom_hetsite_arrays(hetsites:dict):

    chromhetsites = {}
  
    for hetsite in hetsites.values():
        if hetsite.__class__==dict: #dict
            hetsitetype = 'dict'
            chrom = hetsite['chrom']
            if chrom not in chromhetsites:
                chromhetsites[chrom] = []
            chromhetsites[chrom].append(hetsite)
        else: # tuple
            hetsitetype = 'tuple'
            chrom = hetsite.chrom
            if chrom not in chromhetsites:
                chromhetsites[chrom] = []
            chromhetsites[chrom].append(hetsite)
            
    for chrom in chromhetsites:
        if hetsitetype=="dict":
            chromhetsites[chrom].sort(key=lambda h: (h['start'], h['end']))
        else:
            chromhetsites[chrom].sort(key=lambda h: (h.start, h.end))
    
    return chromhetsites

def write_hetallele_bed(hetsitealleles:dict, hetbed:str):

    contigsortedhetalleles = sort_chrom_hetsite_arrays(hetsitealleles)
    logger.debug("Opening " + hetbed + " to write het alleles along assembly contigs")
    with open(hetbed, "w") as hfh:
        for contig in sorted(contigsortedhetalleles.keys()):
            for hetsite in contigsortedhetalleles[contig]:
                hetname = hetsite['name']
                fields = hetname.split("_")
                strand = fields[-1]
                altallele = fields[-2]
                refallele = fields[-3]
                if hetsite['allele'] == refallele:
                    allelehap = 'SAMEHAP'
                elif hetsite['allele'] == altallele:
                    allelehap = 'ALTHAP'
                else:
                    allelehap = 'OTHER'
                assemblycontig = hetsite['query']
                assemblystart = hetsite['start'] - 1
                assemblyend = hetsite['end'] - 1
                hfh.write(contig + "\t" + str(hetsite['start']) + "\t" + str(hetsite['end']) + "\t" + hetsite['name'] + "\t" + hetsite['allele'] + "\t" + hetsite['ref'] + "\t" + str(hetsite['refstart']) + "\t" + str(hetsite['refend']) + "\t" + assemblycontig + "\t" + str(assemblystart) + "\t" + str(assemblyend) + "\t" + allelehap + "\n")

def map_benchmark_hapmers_onto_assembly(queryfasta, matmarkerfile:str, patmarkerfile:str, outputdir:str, outputfiles:dict):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    mathapmeroutput = "KmerMap.mat"
    pathapmeroutput = "KmerMap.pat"
    matpathapmeroutput = "KmerMap.matpat"
    assemblystub = re.sub("\.fa.*", "", queryfasta)
    assemblystub = re.sub(".*/", "", assemblystub)
    matoutputlocation = outputdir + "/" + mathapmeroutput + "." + assemblystub + ".kmers.merge.bed"
    patoutputlocation = outputdir + "/" + pathapmeroutput + "." + assemblystub + ".kmers.merge.bed"
    matpatoutputlocation = outputdir + "/" + matpathapmeroutput + "." + assemblystub + ".hapmers.merge.bed"
    tmpdir = currentdir + "/" + outputdir + "/tmp"
    if not os.path.exists(outputfiles["phasemarkerbed"]):
        if False:
            matpatcommand = "WriteKmerBed -v -m -T2 -P" + tmpdir + " " + matmarkerfile + " " + patmarkerfile + " " + queryfasta + " " + outputdir + "/" + matpathapmeroutput
            path = Path(currentdir + "/" + outputdir + "/tmp")
            logger.info("Creating temporary directory " + outputdir + "/tmp for output")
            path.mkdir(exist_ok=True)
            print("Running: " + matpatcommand)
            logger.info("Running: " + matpatcommand)
            proc = subprocess.Popen(matpatcommand, shell=True, env=env)
            proc.wait()
        else:
            # map maternal/paternal kmers onto the query fasta:
            #matpatoutputlocation = kmers.map_kmer_markers_onto_fasta(queryfasta:str, markerdbs:list, outputdir:str)
            matcommand = "KmerMap -v -m -T2 -P" + tmpdir + " " + matmarkerfile + " " + queryfasta + " " + outputdir + "/" + mathapmeroutput
            patcommand = "KmerMap -v -m -T2 -P" + tmpdir + " " + patmarkerfile + " " + queryfasta + " " + outputdir + "/" + pathapmeroutput
    
            for command in [matcommand, patcommand]:
                path = Path(currentdir + "/" + outputdir + "/tmp")
                logger.info("Creating temporary directory " + outputdir + "/tmp for output")
                path.mkdir(exist_ok=True)
                #print("Running: " + command)
                logger.info("Running: " + command)
                proc = subprocess.Popen(command, shell=True, env=env)
                proc.wait()
    
            logger.debug("Merging output files into " + matpatoutputlocation)
            mergecommand = "cat " + matoutputlocation + " " + patoutputlocation + " > " + matpatoutputlocation
            os.system(mergecommand)
    
        logger.debug("Sorting output file " + matpatoutputlocation)
        sortcommand = "sort -k1,1n -k2,2n -k3,3n " + matpatoutputlocation + " > " + outputfiles["phasemarkerbed"]
        os.system(sortcommand)
    else:
        print("Using output file " + outputfiles["phasemarkerbed"])
        logger.info("Using output file " + outputfiles["phasemarkerbed"])

    [mergedhapmerintervals, mergedhapmerbedfile] = bedtoolslib.mergebed(outputfiles["phasemarkerbed"], collapsecolumns='4', collapseoutput='distinct', collapsedelim=',')

    return 0

def map_benchmark_hapmers_onto_assembly_with_phaseblocks(queryfasta, matmarkerfile:str, patmarkerfile:str, outputdir:str, outputfiles:dict):
    env = os.environ.copy()
    env['LD_LIBRARY_PATH'] = os.getcwd()
    currentdir = os.getcwd()
    phaseblockoutput = "PhaseBlocks"
    assemblystub = re.sub("\.fa.*", "", queryfasta)
    assemblystub = re.sub(".*/", "", assemblystub)
    outputlocation = outputdir + "/" + phaseblockoutput + "." + assemblystub + ".hapmers.bed"
    tmpdir = currentdir + "/" + outputdir + "/tmp"
    command = "PhaseBlocks -v -T2 -P" + tmpdir + " " + matmarkerfile + " " + patmarkerfile + " " + queryfasta + " " + outputdir + "/" + phaseblockoutput
    if not os.path.exists(outputlocation) and not os.path.exists(outputfiles["phasemarkerbed"]):
        #commandargs = command.split(" ")
        path = Path(currentdir + "/" + outputdir + "/tmp")
        logger.info("Creating temporary directory " + outputdir + "/tmp for output")
        path.mkdir(exist_ok=True)
        print("Running: " + command)
        logger.info("Running: " + command)
        proc = subprocess.Popen(command, shell=True, env=env)
        proc.wait()

    if not os.path.exists(outputfiles["phasemarkerbed"]):
        print("Using output file " + outputlocation)
        logger.info("Using output file " + outputlocation)
        os.rename(outputlocation, outputfiles["phasemarkerbed"])
    else:
        print("Using output file " + outputfiles["phasemarkerbed"])
        logger.info("Using output file " + outputfiles["phasemarkerbed"])

    [mergedhapmerintervals, mergedhapmerbedfile] = bedtoolslib.mergebed(outputfiles["phasemarkerbed"], collapsecolumns='4', collapseoutput='distinct', collapsedelim=',')

    return 0

# This routine duplicates the logic of "MerToPhaseBlock" from Arang Rhie's merqury package (but ported from Java to python)
# I've changed some variable names (e.g., "noBreak" to "includegaps") to make it easier (for me) to read:

def find_phase_blocks_from_marker_bed(bedfile:str, scaffnames:list, shortnum=100, shortlimit=20000, includegaps=True):

    [mathap, pathap] = find_haplotype_names_in_hapmer_bedfile(bedfile)

    blockbedstring = ""

    pmathap = re.compile(r'.*MAT.*', re.IGNORECASE)
    matcolor = '255,0,0'
    patcolor = '0,0,255'
    gapcolor = '255,255,255'

    # find scaffold, endpoints, haplotype, and numbers of switches and markers for each phase block
    # also catalog switches (including "Same") in an array to return
    phaseblocks = []
    switches = []

    nummarkers = 0
    numswitches = 0
    numshortswitches = 0
    blockscaff = ""
    blockmarker = ""
    blockstart = -1
    blockend = -1
    shortstart = -1
    shortend = -1
    isshortswitch = True
    with open(bedfile, "r") as bfh:
        markerline = bfh.readline()
        while markerline:
            markerline = markerline.rstrip()
            [scaffold, start, end, markertype] = markerline.split("\t")
            start = int(start)
            end = int(end)

            if markertype == "gap": # gap
                if not includegaps:
                    # create a block with ending at the starting point of the gap, and add it to the list of phaseblocks:
                    blockend = start
                    if blockscaff != "":
                        print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                        #phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                        blockbedstring = blockbedstring + blockscaff + "\t" + str(blockstart) + "\t" + str(blockend) + "\t" + blockmarker + "\t0\t+\t" + str(blockstart) + "\t" + str(blockend) + "\t" + gapcolor + "\n"
                    # start the next block at the end of the gap with "unknown" marker type, without switches and zero markers:
                    blockscaff = scaffnames[int(scaffold)]
                    blockstart = end
                    blockend = end
                    blockmarker = "unknown"
                    numswitches = 0
                    nummarkers = 0
                    isshortswitch = False
            elif (blockscaff == "" or blockscaff != scaffnames[int(scaffold)]): # scaffold has changed or this is start of the first scaff (not a gap line)
                # unless this is the start of the first scaffold, calculate number of markers in the block and add the current block to the list
                if blockscaff != "": # we are currently in a block
                    nummarkers = nummarkers - numshortswitches

                    if pmathap.match(blockmarker):
                        bedcolor = matcolor
                    else:
                        bedcolor = patcolor
                    #print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                    #phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                    blockbedstring = blockbedstring + blockscaff + "\t" + str(blockstart) + "\t" + str(blockend) + "\t" + blockmarker + "\t0\t+\t" + str(blockstart) + "\t" + str(blockend) + "\t" + bedcolor + "\n"
                    # short switches at the ends of scaffolds are included as long range switches
                    if numshortswitches > 0:
                        blockstart = shortstart
                        blockend = shortend
                        nummarkers = numshortswitches
                        blockmarker = markertype
                        if pmathap.match(blockmarker):
                            bedcolor = matcolor
                        else:
                            bedcolor = patcolor
                        print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                        #phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers })
                        blockbedstring = blockbedstring + blockscaff + "\t" + str(blockstart) + "\t" + str(blockend) + "\t" + blockmarker + "\t0\t+\t" + str(blockstart) + "\t" + str(blockend) + "\t" + bedcolor + "\n"
                # start the first block on this new scaffold
                blockmarker = markertype
                numswitches = 0
                numshortswitches = 0
                distfromswitch = 0
                nummarkers = 0
                blockscaff = scaffnames[int(scaffold)]
                blockstart = start
                blockend = end
                if (markertype == mathap or markertype == pathap):
                    nummarkers = nummarkers + 1
                # java version writes out a switch with "Same" here
                #print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Same " + " previous type " + blockmarker)
                #switches.append({'scaffold':scaffold, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Same', 'blockmarker':blockmarker})
                isshortswitch = False
            else: # same scaffold, not a gap
                nummarkers = nummarkers + 1
                # if the marker type matches the block's type (or the last was an unknown at the beginning of a scaffold), extend block end:
                if blockmarker == markertype or blockmarker == "unknown":
                    blockend = end
                    if blockmarker == "unknown":
                        blockmarker = markertype
                    # if we're in a short switch, weve now come back to the block marker type, so reset shortswitch counter, but first add value to number of switches
                    if isshortswitch:
                        numswitches = numswitches + numshortswitches
                        numshortswitches = 0
                        isshortswitch = False
                        # java version writes out a switch with "SwitchBack" here
                        #print("SWITCH " + blockscaff + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch SwitchBack " + " previous type " + blockmarker)
                        #switches.append({'scaffold':blockscaff, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'SwitchBack', 'blockmarker':blockmarker})
                    #else:
                        # java version writes out a switch with "Same" here
                        #print("SWITCH " + scaffold + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Same " + " previous type " + blockmarker)
                        #switches.append({'scaffold':blockscaff, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Same', 'blockmarker':blockmarker})
                # if the marker is the opposite of the block's type, extend the short switch and increase the counter
                elif blockmarker != markertype:
                    shortend = end
                    numshortswitches = numshortswitches + 1
                    # if we were not in a short switch already, set the start of the switch region to this marker's start:
                    if not isshortswitch:
                        shortstart = start
                        isshortswitch = True
                        # java version writes out a switch with "Short" here
                        #print("SWITCH " + blockscaff + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Short " + " previous type " + blockmarker)
                        #switches.append({'scaffold':blockscaff, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Short', 'blockmarker':blockmarker})
                    # if we were already in a short switch, check to see if it's numerous or long enough to be a long switch:
                    else:
                        distfromswitch = shortend - shortstart
                        # if we're in a long switch, decrease the number of markers by the number of short switches, and add a block
                        if numshortswitches >= shortnum or distfromswitch > shortlimit:
                            nummarkers = nummarkers - numshortswitches
                            print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(blockend) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
                            #phaseblocks.append({'scaffold':blockscaff, 'blockstart':blockstart, 'blockend':blockend, 'blockmarker':blockmarker, 'numswitches':numswitches, 'nummarkers':nummarkers})
                            if pmathap.match(blockmarker):
                                bedcolor = matcolor
                            else:
                                bedcolor = patcolor
                            blockbedstring = blockbedstring + blockscaff + "\t" + str(blockstart) + "\t" + str(blockend) + "\t" + blockmarker + "\t0\t+\t" + str(blockstart) + "\t" + str(blockend) + "\t" + bedcolor + "\n"
                            # then initiate a new block on the switched haplotype (allowing the formerly short block to extend)
                            blockstart = shortstart
                            blockend = end
                            nummarkers = numshortswitches
                            numswitches = 0
                            numshortswitches = 0
                            blockmarker = markertype
                            isshortswitch = False
                            # java version writes out a switch with "Long" here
                            #print("SWITCH " + blockscaff + ":" + str(start) + "-" + str(end) + " type " + markertype + " switch Long " + " previous type " + blockmarker)
                            #switches.append({'scaffold':blockscaff, 'start':start, 'end':end, 'markertype':markertype, 'switchtype':'Long', 'blockmarker':blockmarker})
            markerline = bfh.readline()

        #print("Appending block " + blockscaff + ":" + str(blockstart) + "-" + str(end) + " type " + blockmarker + " with " + str(numswitches) + " switches and " + str(nummarkers) + " markers")
        if pmathap.match(blockmarker):
            bedcolor = matcolor
        else:
            bedcolor = patcolor
        blockbedstring = blockbedstring + blockscaff + "\t" + str(blockstart) + "\t" + str(end) + "\t" + blockmarker + "\t0\t+\t" + str(blockstart) + "\t" + str(end) + "\t" + bedcolor + "\n"

    phaseblockbed = pybedtools.BedTool(blockbedstring, from_string = True)

    return phaseblockbed

# Use the Viterbi algorithm to determine the most probable chain of phase block haplotype states #
# given a set of observed haplotype-specific kmers from the benchmark along the test assemblies' scaffolds #
# scaffold fasta object should be passed if you want phase block coordinates to extend to the ends of the 
# scaffolds--otherwise they will only extend to the outermost marker positions!
def find_hapmer_phase_blocks_with_hmm(bedfile:str, hmmphaseblockbedfile:str, scafffastaobj:str, alpha:float, beta:float, distmultiplier=0 ):

    emissionprobs = [[1.0-alpha, alpha], [alpha, 1.0-alpha]]

    # determine maternal/paternal haplotype names in marker bed--these could be mat/pat or hap1/hap2
    [mathap, pathap] = find_haplotype_names_in_hapmer_bedfile(bedfile)

    if mathap is None:
        logger.info("Only haplotype found in " + bedfile + " is " + pathap + "--returning entirely one state")
        mathap = "missinghap"
    if pathap is None:
        logger.info("Only haplotype found in " + bedfile + " is " + mathap + "--returning entirely one state")
        pathap = "missinghap"

    if mathap == "missinghap" and pathap == "missinghap":
        logger.critical("No haplotype markers found for phasing chromosomes!")
        exit(1)

    shortmathap = re.sub("_{0,1}not.*", "", mathap)
    shortmathap = re.sub(".*\.", "", shortmathap)
    shortpathap = re.sub("_{0,1}not.*", "", pathap)
    shortpathap = re.sub(".*\.", "", shortpathap)
    logger.debug("Using " + shortmathap + " for first haplotype name and " + shortpathap + " for second")
    phaseblockbedstring = ""
    mathap1color = '255,0,0'
    pathap2color = '0,0,255'

    haplogprobs = []
    # time 0 has no previous state values:
    previousstate = [[None, None]]
    timepoint = 0

    current_scaffold = ""
    scaffname = ""
    current_pos = 0
    intermarkercoords = []
    phasedend1 = "." + mathap + ".bed"
    phasedend2 = "." + pathap + ".bed"
    matphaseblockbed = hmmphaseblockbedfile.replace('.bed', phasedend1)
    patphaseblockbed = hmmphaseblockbedfile.replace('.bed', phasedend2)
    scaffnamemarkerbed = bedfile.replace('.bed', '.scaffnames.bed')
    mfh = open(matphaseblockbed, "w")
    pfh = open(patphaseblockbed, "w")
    sfh = open(scaffnamemarkerbed, "w")
    scaffnames = scafffastaobj.references
    with open(bedfile, "r") as bfh:
        markerline = bfh.readline()
        while markerline:
            markerline = markerline.rstrip()
            [scaffold, start, end, markertype] = markerline.split("\t")
            scaffold = int(scaffold)
            start = int(start)
            end = int(end)
            midpos = int((start + end)/2)
            scaffname = scaffnames[scaffold]

            # if new midpoint is smaller than current position and on same scaffold, skip it
            # (this is usually because of overlapping regions from different haplogroups)
            if current_scaffold == scaffold and midpos <= current_pos:
                markerline = bfh.readline()
                continue

            if markertype == pathap:
                bedcolor = pathap2color
            else:
                bedcolor = mathap1color
            sfh.write(scaffname + "\t" + str(start) + "\t" + str(end) + "\t" + markertype + "\t0\t+\t" + str(start) + "\t" + str(end) + "\t" + bedcolor + "\n")

            if current_scaffold != scaffold:
                # write out max log prob states if there was a previous chromosome:
                if current_scaffold != "":
                    if len(haplogprobs) > 0:
                        oldscaffname = scaffnames[current_scaffold]
                        lastlogprobs = haplogprobs.pop()
                        scafflength = scafffastaobj.get_reference_length(oldscaffname)
                        write_scaffold_phase_blocks(oldscaffname, scafflength, timepoint, lastlogprobs, previousstate, intermarkercoords, mathap1color, pathap2color, shortmathap, shortpathap, mfh, pfh)
                    timepoint = 0
                current_scaffold = scaffold
                current_pos = 0
                intermarkercoords = []
                previousstate = [[None, None]]
                haplogprobs = []
            if markertype == mathap or markertype == pathap:
                if markertype == mathap:
                    hapindex = 0
                else:
                    hapindex = 1
                # for each timepoint, append an entry to the array of state probabilities--it will be overwritten
                haplogprobs.append([None, None])
                if timepoint == 0:
                    haplogprobs[0][0] = math.log(emissionprobs[0][hapindex])
                    haplogprobs[0][1] = math.log(emissionprobs[1][hapindex])
                    #logger.debug("Assigning haploprobs " + str(haplogprobs[timepoint][0]) + " and " + str(haplogprobs[timepoint][1]) + " at timepoint " + str(timepoint))
                else:
                    previousstate.append([0,0])
                    transitionprobs = [[1.0-beta, beta], [beta, 1.0-beta]]
                    for state in range(2):
                        for laststate in range(2):
                            transitionprob = transitionprobs[laststate][state]
                            newlogprob = haplogprobs[timepoint-1][laststate] + math.log(transitionprob) + math.log(emissionprobs[state][hapindex])
                            if haplogprobs[timepoint][state] is None or newlogprob > haplogprobs[timepoint][state]:
                                haplogprobs[timepoint][state] = newlogprob
                                previousstate[timepoint][state] = laststate
                                #print("Logging previous state " + str(laststate) + " at time " + str(timepoint))
                    #logger.debug("Assigning haploprobs " + str(haplogprobs[timepoint][0]) + " and " + str(haplogprobs[timepoint][1]) + " at timepoint " + str(timepoint))
                intermarkercoords.append([current_pos, midpos])
                current_pos = midpos
                timepoint = timepoint + 1
            markerline = bfh.readline()

    scafflength = scafffastaobj.get_reference_length(scaffname)
    write_scaffold_phase_blocks(scaffname, scafflength, timepoint, lastlogprobs, previousstate, intermarkercoords, mathap1color, pathap2color, shortmathap, shortpathap, mfh, pfh)

    mfh.close()
    pfh.close()
    sfh.close()

    [matblockintervals, matblockbedfile] = bedtoolslib.mergebed(matphaseblockbed, collapsecolumns='4,5,6,7,8,9', collapseoutput='distinct,first,first,first,last,first', collapsedelim=',')
    [patblockintervals, patblockbedfile] = bedtoolslib.mergebed(patphaseblockbed, collapsecolumns='4,5,6,7,8,9', collapseoutput='distinct,first,first,first,last,first', collapsedelim=',')

    allphaseblockints = matblockintervals.cat(patblockintervals, postmerge=False).sort()

    return allphaseblockints

def find_haplotype_names_in_hapmer_bedfile(bedfile:str):

    # if one found haplotype matches ".mat" or "[ab]1_not_[ab]2" (case insensitive), return it first--covers both benchmark and assembly comparison
    pmathap = re.compile(r'.*\.mat.*|[qr]1_not_[qr]2', re.IGNORECASE)
    pcomma = re.compile(r'.*\,.*', re.IGNORECASE)
    # determine two haplotypes in the "name" field of the input bed file (assumes only two distinct values plus optionally "gap")
    with open(bedfile, "r") as bfh:
        markerline = bfh.readline()
        hap1 = None
        hap2 = None
        while markerline:
            markerline = markerline.rstrip()
            [scaffold, start, end, markertype] = markerline.split("\t")
            if pcomma.match(markertype):
                logger.debug("Skipping marker name " + markertype)
                markerline = bfh.readline()
                continue
            if markertype != "gap" and hap1 is None:
                hap1 = markertype
                logger.debug("Found first marker name " + markertype)
            elif markertype != "gap" and markertype != hap1:
                hap2 = markertype
                logger.debug("Found second marker name " + markertype)
                break
            markerline = bfh.readline()
    if hap1 is not None and hap2 is not None:
        if pmathap.match(hap1):
            logger.info("Found hap1  " + hap1 + " and hap2 " + hap2 + " as haplotypes.")
            print("Found hap1 " + hap1 + " and hap2 " + hap2 + " as haplotypes.")
            return [hap1, hap2]
        elif pmathap.match(hap2):
            logger.info("Found hap2 " + hap2 + " and hap1  " + hap1 + " as haplotypes.")
            print("Found hap2 " + hap2 + " and  hap1 " + hap1 + " as haplotypes.")
            return [hap2, hap1]
        else:
            logger.info("Unable to determine mat/pat or 1/2 identities of hap1 " + hap1 + " and hap2 " + hap2)
            print("Unable to determine mat/pat or 1/2 identities of hap1 " + hap1 + " and hap2 " + hap2)
            exit(1)
    else:
        logger.info("Unable to find two haplotypes in " + bedfile)
        print("Unable to find two haplotypes in " + bedfile)
        if hap1 is not None:
            if pmathap.match(hap1):
                return[hap1, hap2]
            else:
                return [hap2, hap1]
        if hap2 is not None:
            if pmathap.match(hap2):
                return [hap2, hap1]
            else:
                return [hap1, hap2]

        return [hap1, hap2]

def write_scaffold_phase_blocks(scaffoldname:str, scafflength:int, timepoint:int, lastlogprobs:list, previousstate:list, intermarkercoords:list, hap1color:str, hap2color:str, mathap:str, pathap:str, hap1fh, hap2fh):
    prevstatelength = len(previousstate)
    intermarkerlength = len(intermarkercoords)
    path = [0] * timepoint
    if lastlogprobs[0] > lastlogprobs[1]:
        path[timepoint-1] = 0
    else:
        path[timepoint-1] = 1
    revtimepoints = range(timepoint-2, -1, -1)
    for currenttime in revtimepoints:
        path[currenttime] = previousstate[currenttime + 1][path[currenttime + 1]]
    for currenttime in range(timepoint):
        if path[currenttime] == 0:
            currenthap = mathap
            fh = hap1fh
            bedcolor = hap1color
        else:
            currenthap = pathap
            fh = hap2fh
            bedcolor = hap2color

        if currenttime == 0:
            blockstart = 0
        else:
            blockstart = int((intermarkercoords[currenttime][0] + intermarkercoords[currenttime][1])/2.0)

        if currenttime < timepoint-1:
            blockend = int((intermarkercoords[currenttime+1][0] + intermarkercoords[currenttime+1][1])/2.0)
        else:
            blockend = scafflength
        fh.write(scaffoldname + "\t" + str(blockstart) + "\t" + str(blockend) + "\t" + currenthap + "\t0\t+\t" + str(blockstart) + "\t" + str(blockend) + "\t" + bedcolor + "\n")

    return 0
