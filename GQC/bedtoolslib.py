import os
import re
import pybedtools
import pysam
import logging

logger = logging.getLogger(__name__)

def mergebed(bedfile:str, collapsecolumns='4', collapseoutput='collapse', collapsedelim='|', skipifexists=True)->str:

    mergedbed = bedfile.replace(".bed", ".merged.bed")

    if not skipifexists or not os.path.exists(mergedbed):
        logger.debug("Merging " + bedfile + " to create " + mergedbed)

        unmergedints = pybedtools.bedtool.BedTool(bedfile)
        numunmerged = unmergedints.count()
   
        logger.debug("There are " + str(numunmerged) + " unmerged intervals")

        if numunmerged > 0:
            mergedints = unmergedints.merge(c=collapsecolumns, o=collapseoutput, delim=collapsedelim)
            mergedints.saveas(mergedbed)
        else:
            with open(mergedbed, 'a'):
                os.utime(mergedbed, None) 
            mergedints = unmergedints
    else:
        logger.debug("Skipping merging of " + bedfile + ": merged bedfile " + mergedbed + " already exists")
        mergedints = pybedtools.bedtool.BedTool(mergedbed)

    return [mergedints, mergedbed]

def mergeintervals(intervals, collapsecolumns='4', collapseoutput='collapse', collapsedelim='|', distance=0):

    mergedints = intervals.merge(c=collapsecolumns, o=collapseoutput, delim=collapsedelim, d=distance)

    return mergedints

def mergemultiplebedfiles(bedfilelist:list, sort=True, postmerge=False):

    if len(bedfilelist) < 2:
        logger.critical("Cannot call mergemultiplebedfiles on less than two bed files!")
        exit(1)
    if len(bedfilelist) > 2:
        logger.critical("Theres a bug in GQC bedtoolslib.py so it cant yet merge multiple excluded regions bed files!")
        exit(1)

    firstbedtool = pybedtools.bedtool.BedTool(bedfilelist[0])
    #allbedtools = firstbedtool.cat(bedfilelist[1:], postmerge=postmerge)
    allbedtools = firstbedtool.cat(bedfilelist[1])

    if sort:
        return allbedtools.sort()
    else:
        return allbedtools

def bedsum(intervals)->int:

    alllengths = map(len, intervals)
    bedsum = sum(alllengths)

    return bedsum

def intersectbed(bedfile1:str, bedfile2:str, outputfile:str, writefirst=False, writeboth=False, outerjoin=False)->str:

    logger.debug("Intersecting " + bedfile1 + " and " + bedfile2 + " to create " + outputfile)

    command = "bedtools intersect -a " + bedfile1 + " -b " + bedfile2
    if writefirst:
        command = command + " -wa"
    if writeboth:
        command = command + " -wo"
    if outerjoin:
        command = command + " -loj"
    os.system(command + " > " + outputfile)
    intersectbed = pybedtools.bedtool.BedTool(outputfile)

    return [intersectbed, outputfile]

def intersectintervals(intervals1, intervals2, v=False, wo=False, wa=False, wb=False, counts=False):

    if wo:
        print("Intersect called with -wo option")
    intersectedints = intervals1.intersect(intervals2, v=v, wa=wa, wb=wb, wo=wo, c=counts)

    return intersectedints

def subtractintervals(intervals1, intervals2, A=False, wb=False, wo=False):

    subtractedints = intervals1.subtract(intervals2, A=A, wb=wb, wo=wo)

    return subtractedints

def genomeintervals(fastafile:str):

    indexfile = fastafile + ".fai"
    if not os.path.exists(indexfile):
        pysam.faidx(fastafile)

    fastaobj = pysam.FastaFile(fastafile)
    bedstring = ""
    for refname in fastaobj.references:
        reflength = fastaobj.get_reference_length(refname)
        bedstring = bedstring + refname + "\t0\t" + str(reflength) + "\n"
    
    return pybedtools.BedTool(bedstring, from_string=True)

