import os
import re
import sys
import random
import pybedtools
import logging
import datetime
from collections import namedtuple
from GQC import seqparse
from GQC import alignparse
from GQC import bedtoolslib
from GQC import coverage

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

logger = logging.getLogger(__name__)

def classify_errors(refobj, queryobj, variants, hetsites, outputdict, benchparams, stats, args)->str:

    bencherrorfile = outputdict["bencherrortypebed"]
    stats["singlebasecounts"] = {}
    stats["indellengthcounts"] = {}
    stats["totalerrorsinaligns"] = 0

    testerrorfile = outputdict["testerrortypebed"]
    tfh = open(testerrorfile, "w")

    xfh = None
    if "benchexcludederrortypebed" in outputdict.keys():
        excludederrorfile = outputdict["benchexcludederrortypebed"]
        xfh = open(excludederrorfile, "w")
    if args.vcf:
        bencherrorvcf = bencherrorfile.replace(".bed", "")
        bencherrorvcf = bencherrorvcf + ".vcf"
        vfh = open(bencherrorvcf, "w")
        vcfheader = vcf_header(args)
        vfh.write(vcfheader)
        if "assembly" in args:
            samplename = args.assembly
        elif "qname" in args:
            samplename = args.qname
        else:
            samplename = "SAMPLE"

        vfh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + samplename + "\n")

    with open(bencherrorfile, "w") as efh:
        for variant in variants:
            namefields = variant.name.split("_")
            numfields = len(namefields)
            contigname = "_".join(namefields[0:numfields-4])
            pos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]
            if namefields[-1] == "F":
                alignstrand = '+'
            else:
                alignstrand = '-'
            #logger.debug("Splitting variant name " + variant.name + " and found contigname " + contigname)
            varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele

            if hetsites and varname in hetsites.keys():
                errortype = 'PHASING'
                errortypecolor = '255,0,0'
            else:
                errortype = 'CONSENSUS'
                errortypecolor = '0,0,255'

            benchvarstart = variant.start
            benchvarend = variant.end
            if refallele == "*":
                benchvarstart = benchvarstart - 1

            varqvscore = "1000"
            if variant.qvscore is not None:
                varqvscore = str(variant.qvscore)

            if xfh and variant.excluded:
                xfh.write(variant.chrom + "\t" + str(benchvarstart) + "\t" + str(benchvarend) + "\t" + varname + "\t" + varqvscore + "\t" + alignstrand + "\t" + str(benchvarstart) + "\t" + str(benchvarend) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + variant.name + "\n")
                if args.vcf:
                    vcfrecord = vcf_format(variant, refobj, queryobj)
                    vfh.write(vcfrecord)
            else:
                efh.write(variant.chrom + "\t" + str(benchvarstart) + "\t" + str(benchvarend) + "\t" + varname + "\t" + varqvscore + "\t" + alignstrand + "\t" + str(benchvarstart) + "\t" + str(benchvarend) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + variant.name + "\n")
                if args.vcf:
                    vcfrecord = vcf_format(variant, refobj, queryobj)
                    vfh.write(vcfrecord)

            tfh.write(contigname + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + variant.name + "\t" + varqvscore + "\t" + alignstrand + "\t" + str(pos-1) + "\t" + str(pos - 1 + len(altallele)) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.vartype + "\t" + varname + "\n")
            stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                
            # tally statistics for non-phasing errors:

            if variant.excluded or errortype != "CONSENSUS":
                continue
            if variant.vartype == "SNV":
                snvkey = refallele + "_" + altallele
                if snvkey in stats["singlebasecounts"]:
                    stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                else:
                    stats["singlebasecounts"][snvkey] = 1
            else:
                trueref = refallele.replace('*', '')
                truealt = altallele.replace('*', '')
                lengthdiff = len(truealt) - len(trueref)
                if lengthdiff in stats["indellengthcounts"]:
                    stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                else:
                    stats["indellengthcounts"][lengthdiff] = 1
    tfh.close()
    if xfh:
        xfh.close()

    if args.vcf:
        vfh.close()

    return stats

def vcf_format(variant, refobj, queryobj):

    chrom = variant.chrom
    start = variant.start # zero-based
    refpos = start + 1
    namefields = variant.name.split("_")
    numfields = len(namefields)
    contigname = "_".join(namefields[0:numfields-4])
    contigpos = int(namefields[-4])
    refallele = namefields[-3]
    altallele = namefields[-2]
    if namefields[-1] == "F":
        alignstrand = '+'
    else:
        alignstrand = '-'
    filterfield = 'PASS'
    if variant.excluded:
        filterfield = 'EXCLUDED'

    refallele = refallele.replace("*", "")
    altallele = altallele.replace("*", "")
    while len(refallele) > 0 and len(altallele) > 0 and refallele[-1] == altallele[-1]:
        refallele = refallele[:-1]
        altallele = altallele[:-1]

    if refallele == "" or altallele == "":
        refpos = refpos - 1
        refallele = refobj.fetch(reference=chrom, start=refpos-1, end=refpos) + refallele
        refallele = refallele.upper()
        if alignstrand == '+':
            contigpos = contigpos - 1
            altallele = queryobj.fetch(reference=contigname, start=contigpos-1, end=contigpos) + altallele
            altallele = altallele.upper()
        else:
            if len(refallele) == 1:
                contigend = contigpos + len(altallele)
            else:
                contigend = contigpos + len(altallele) + 1
            altbase = queryobj.fetch(reference=contigname, start=contigend-1, end=contigend)
            altbase = seqparse.revcomp(altbase)
            altbase = altbase.upper()
            altallele = altbase + altallele
            logger.debug("Reverse strand empty allele adjustment: fetch query " + contigname + ":" + str(contigend) + "-" + str(contigend) + " gives " + refallele + "/" +  altallele)
            logger.debug(chrom + "\t" + str(refpos) + "\t" + variant.name + "\t" + refallele + "\t" + altallele + "\t.\t"  + filterfield + "\t.\tGT\t1")

    return chrom + "\t" + str(refpos) + "\t" + variant.name + "\t" + refallele + "\t" + altallele + "\t.\t"  + filterfield + "\t.\tGT\t1\n"

def vcf_header(args):

    if 'benchmark' in args:
        benchname = args.benchmark
    elif 'rname' in args:
        benchname = args.rname
    else:
        benchname = 'REF'
    date = datetime.datetime.now()
    datestring = date.strftime("%Y%m%d")
    header_string = "##fileformat=VCFv4.5\n##fileDate=" + datestring + "\n##source=GQC\n##reference=" + benchname + "\n##FILTER=<ID=EXCLUDED, Description=\"In excluded region of the benchmark reference\">\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"

    return header_string

def gather_mononuc_stats(coveredmononucbedfile:str, mononucstatsfile:str):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")
    result = {}
    sfh = open(mononucstatsfile, "w")
    with open(coveredmononucbedfile, "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color, chrom_b, start_b, end_b, variant_name, score_b, strand_b, widestart_b, wideend_b, color_b, error_type, variant_type, queryvariantname] = mononucline.split("\t")
            runlength = int(end) - int(start)
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            if chrom_b == '.':
                result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':runlength, 'type':'CORRECT'}
                sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(runlength) + "\tCORRECT\n")
            else:
                variantnamefields = variant_name.split("_")
                refbases = variantnamefields[-2]
                altbases = variantnamefields[-1]
                if refbases == "*" and p[repeatedbase].match(altbases): # increased length
                    newlength = runlength + len(altbases)
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                elif altbases == "*" and p[repeatedbase].match(refbases): # decreased length
                    newlength = runlength - len(refbases)
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                elif p[repeatedbase].match(refbases) and p[repeatedbase].match(altbases): # expanded notation
                    newlength = len(altbases)
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                else: # complex error
                    newlength = -1
                    result[name] = {'base':repeatedbase, 'length':runlength, 'assemblylength':newlength, 'type':error_type}
                    sfh.write(name + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(newlength) + "\t" + error_type + "\n")

            mononucline = mfh.readline()

    return result

def assess_mononuc_read_coverage(align_obj, mononucbedfile, outputdict, bedintervals, hetsitedict, args):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")

    max_reads_per_mononuc = args.maxmononucreads

    mononucstatsfile = outputdict["mononucstatsfile"]

    mononucbedintervals = pybedtools.BedTool(mononucbedfile)
    if bedintervals is not None:
        includedmononucbeds = bedtoolslib.intersectintervals(mononucbedintervals, bedintervals, wa=True)
        includedmononucbeds.saveas(outputdict["includedmononucfile"])
    else:
        mononucbedintervals.saveas(outputdict["includedmononucfile"])

    mononucdict = {}
    msfh = open(mononucstatsfile, "w", buffering=1)
    # for each mononuc run in the benchmark, retrieve reads, count length of mononuc, classify as CORRECT, HET, ERROR, or COMPLEX
    with open(outputdict["includedmononucfile"], "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color] = mononucline.split("\t")
            runlength = int(end) - int(start)
            if runlength not in mononucdict:
                mononucdict[runlength] = {}
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            refalleleseq = repeatedbase * runlength
            readsforthismononuc = 0
            for readalign in align_obj.fetch(contig=chrom, start=int(start), stop=int(end)):
                if max_reads_per_mononuc > 0 and readsforthismononuc >= max_reads_per_mononuc:
                    break
                if args.downsample is not None and random.random() >= args.downsample:
                    continue
                queryseq = readalign.query_sequence
                if readalign.is_secondary or queryseq is None:
                    continue
                readname = readalign.query_name
                if readalign.reference_start >= int(start) or readalign.reference_end <= int(end):
                    continue
                pairs = readalign.get_aligned_pairs()
                if len(pairs) == 0:
                    continue
                # find zero-based read pos of base aligned to ref base one before mononuc:
                readstart = find_readpos_in_pairs(pairs, int(start)-1)
                # find zero-based read pos of base aligned to ref base one after mononuc:
                readend = find_readpos_in_pairs(pairs, int(end))
                if readstart is not None and readend is not None:
                    queryseq = queryseq.upper()
                    alleleseq = queryseq[readstart+1:readend]
                    # extend the read sequence so 5' and 3' to include all of repeated bases
                    #if readstart >=0:
                        #logger.debug("Checking " + queryseq[readstart] + " vs " + repeatedbase)
                    while readstart >= 0 and queryseq[readstart] == repeatedbase:
                        #logger.debug("Extending " + alleleseq + " one base")
                        alleleseq = repeatedbase + alleleseq
                        readstart = readstart - 1
                    #if readend < len(queryseq):
                        #logger.debug("Checking " + queryseq[readend] + " vs " + repeatedbase)
                    while readend < len(queryseq) and queryseq[readend] == repeatedbase:
                        #logger.debug("Extending " + alleleseq + " one base")
                        alleleseq = alleleseq + repeatedbase
                        readend = readend + 1

                    if p[repeatedbase].match(alleleseq):
                        numbases = len(alleleseq)
                        matchtype = "CORRECT"
                        if numbases != runlength:
                            potentialhetname = chrom + "_" + str(int(start)+1) + "_" + refalleleseq + "_" + alleleseq
                            if potentialhetname in hetsitedict:
                                matchtype = "HET"
                            else:
                                matchtype = "ERROR"
                    else:
                        numbases = -1
                        matchtype = "COMPLEX"

                    msfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(numbases) + "\t" + matchtype + "\n")
                    if numbases not in mononucdict[runlength]:
                        mononucdict[runlength][numbases] = {}
                    if matchtype not in mononucdict[runlength][numbases]:
                        mononucdict[runlength][numbases][matchtype] = 1
                    else:
                        mononucdict[runlength][numbases][matchtype] = mononucdict[runlength][numbases][matchtype] + 1
                    readsforthismononuc = readsforthismononuc + 1
                else:
                    logger.debug(readname + " is unaligned at one endpoint of interval " + str(start) + " to " + str(end) + "!")
                    if readstart is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a start")
                    if readend is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a end")

            mononucline = mfh.readline()

    return mononucdict

def assess_str_read_coverage_require_flank(strlabel, align_obj, refobj, strbedfile, outputdict, bedintervals, hetsitedict, args):

    max_reads_per_str = args.maxstrreads

    includedstrfile = outputdict["includedstrprefix"] + "." + strlabel + ".included.bed"
    strstatsfile = outputdict["strstatsprefix"] + "." + strlabel + ".stats.out"

    logger.info("Writing " + strlabel + " STR data to " + strstatsfile)

    # if the file already exists, don't recreate it
    if os.path.exists(strstatsfile) and os.path.getsize(strstatsfile) > 0:
        strdict = retrieve_stats_from_strstats_file(strstatsfile)
        return strdict

    strbedintervals = pybedtools.BedTool(strbedfile)
    if bedintervals is not None:
        includedstrbeds = bedtoolslib.intersectintervals(strbedintervals, bedintervals, wa=True)
        includedstrbeds.saveas(includedstrfile)
    else:
        strbedintervals.saveas(includedstrfile)

    print("Saved included STR bed file. Will now open statsfile")
    strdict = {}
    sfh = open(strstatsfile, "w", buffering=1)
    # for each str run in the str bedfile, retrieve reads, count length of str, classify as CORRECT, HET, ERROR, or COMPLEX
    with open(includedstrfile, "r") as ifh:
        strline = ifh.readline()
        while strline:
            strline = strline.rstrip()
            strfields = strline.split("\t")
            chrom = strfields[0]
            start = int(strfields[1])
            end = int(strfields[2])
            name = strfields[3]

            # need to extract benchmark's flanking bases here for comparison later!!!
            benchleftflankseq = refobj.fetch(reference=chrom, start=start-5, end=start)
            benchleftflankseq = benchleftflankseq.upper()
            benchrightflankseq = refobj.fetch(reference=chrom, start=end, end=end+5)
            benchrightflankseq = benchrightflankseq.upper()
            runlength = end - start
            if runlength not in strdict:
                strdict[runlength] = {}
            namefields = name.split("_")
            repeatedbases = namefields[-1]
            refalleleseq = refobj.fetch(reference=chrom, start=start, end=end)
            readsforthisstr = 0
            # for each read, attempt to find read positions aligned to the bases immediately flanking
            # the STR. If one side and/or the other has no aligning base (due either to a deletion or
            # to clipped sequence) then check for clipping beginning at that base, and
            # extend the STR-aligned read sequence with clipped sequence to include adjacent 
            # bases that extend the repeat plus 5 bp additional (this is NIST's so-called "slop" length).
            # If a flanking base *does* have a read base aligned to it, also extend that end to 
            # include repeated bases plut 5 bp slop.
            for readalign in align_obj.fetch(contig=chrom, start=start, stop=end):
                if max_reads_per_str > 0 and readsforthisstr >= max_reads_per_str:
                    break
                if args.downsample is not None and random.random() >= args.downsample:
                    continue
                # query_sequence includes soft clipped bases:
                queryseq = readalign.query_sequence
                if readalign.is_secondary or queryseq is None:
                    continue
                queryseq = queryseq.upper()
                readname = readalign.query_name
                # if the read doesn't align to the entire STR, skip it:
                if readalign.reference_start is None or readalign.reference_end is None:
                    logger.debug("Read " + readname + " has no reference_start or reference_end for " + chrom + ":" + str(start) + "-" + str(end))
                    continue
                if readalign.reference_start > start or readalign.reference_end < end:
                    continue
                # if the read doesn't align to the flanking bases, it must be clipped:
                cigartuples = readalign.cigartuples
                if ((readalign.reference_start == start and cigartuples[0][0] != 4) or
                    (readalign.reference_end == end and cigartuples[len(cigartuples)-1][0] != 4)):
                    logger.debug("Read " + readname + " aligns just to end of " + name + " but is not soft clipped--skipping")
                    continue
                pairs = readalign.get_aligned_pairs()
                if len(pairs) == 0:
                    continue
                # find zero-based read pos of base aligned to ref base one before STR:
                if readalign.reference_start < start:
                    readstart = find_readpos_in_pairs(pairs, start-1)
                elif readalign.query_alignment_start > 0: # must be clipped--set start to one efore base aligning to reference start:
                    readstart = readalign.query_alignment_start - 1
                # find zero-based read pos of base aligned to ref base one after STR:
                if readalign.reference_end > end:
                    readend = find_readpos_in_pairs(pairs, end)
                elif readalign.query_alignment_end < len(queryseq):
                    readend = readalign.query_alignment_end + 1

                if readstart is not None and readend is not None:
                    # extend the read sequence to 5' and 3' so long as it matches the STR. "correctcomp" is True or False depending
                    # on whether the alleleseq is a faithful repeat (even if flanking copies are partial) of the repeatedbases passed to it:
                    [alleleseq, readstart, readend, correctcomp] = widen_str_on_readseq(queryseq, readstart, readend, repeatedbases)

                    # check that read has 5 flanking base pairs on either side of the HP:
                    if readstart < 5 or readend > len(queryseq)-5:
                        logger.debug("Skipping align of read " + readname + " which has fewer than 5 bp outside of HP repeat for " + chrom + ":" + str(start) + "-" + str(end))
                        continue
                    else:
                        leftflankseq = queryseq[readstart-4:readstart+1]
                        leftflankseq = leftflankseq.upper()
                        rightflankseq = queryseq[readend:readend+5]
                        rightflankseq = rightflankseq.upper()
                    matchtype = "UNKNOWN"
                    if correctcomp:
                        numbases = len(alleleseq)
                        if numbases != runlength:
                            potentialhetname = chrom + "_" + str(start+1) + "_" + refalleleseq + "_" + alleleseq
                            if potentialhetname in hetsitedict:
                                matchtype = "HET"
                            else:
                                matchtype = "LENGTHERROR"
                                logger.debug("Potential het " + potentialhetname + " is not a het")
                        else:
                            if leftflankseq == benchleftflankseq and rightflankseq == benchrightflankseq:
                                matchtype = "CORRECT"
                            else:
                                logger.debug(leftflankseq + "/" + rightflankseq + " doesnt match " + benchleftflankseq + "/" + benchrightflankseq)
                                matchtype = "FLANKERROR"
                    else: # discrepancies within the HP 
                        numbases = -1
                        matchtype = "COMPLEX"

                    sfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbases + "\t" + str(runlength) + "\t" + str(numbases) + "\t" + matchtype + "\n")
                    if numbases not in strdict[runlength]:
                        strdict[runlength][numbases] = {}
                    if matchtype not in strdict[runlength][numbases]:
                        strdict[runlength][numbases][matchtype] = 1
                    else:
                        strdict[runlength][numbases][matchtype] = strdict[runlength][numbases][matchtype] + 1
                    readsforthisstr = readsforthisstr + 1
                else:
                    logger.debug(readname + " is unaligned at one endpoint of interval " + str(start) + " to " + str(end) + " and has fewer than 5 clipped bases!")
                    if readstart is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a start")
                    if readend is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a end")

            strline = ifh.readline()

    return strdict

def assess_mononuc_read_coverage_require_flank(align_obj, refobj, mononucbedfile, outputdict, bedintervals, hetsitedict, args):

    p = {}
    p["A"] = re.compile("^[aA]+$")
    p["T"] = re.compile("^[tT]+$")
    p["C"] = re.compile("^[cC]+$")
    p["G"] = re.compile("^[gG]+$")

    max_reads_per_mononuc = args.maxmononucreads

    mononucstatsfile = outputdict["mononucstatsfile"]

    # if the file already exists, don't recreate it
    if os.path.exists(mononucstatsfile):
        mononucdict = retrieve_stats_from_strstats_file(mononucstatsfile)
        return mononucdict

    mononucbedintervals = pybedtools.BedTool(mononucbedfile)
    if bedintervals is not None:
        includedmononucbeds = bedtoolslib.intersectintervals(mononucbedintervals, bedintervals, wa=True)
        includedmononucbeds.saveas(outputdict["includedmononucfile"])
    else:
        mononucbedintervals.saveas(outputdict["includedmononucfile"])

    mononucdict = {}
    msfh = open(mononucstatsfile, "w", buffering=1)
    # for each mononuc run in the benchmark, retrieve reads, count length of mononuc, classify as CORRECT, HET, ERROR, or COMPLEX
    with open(outputdict["includedmononucfile"], "r") as mfh:
        mononucline = mfh.readline()
        while mononucline:
            mononucline = mononucline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color] = mononucline.split("\t")
            # need to extract benchmark's flanking bases here for comparison later!!!
            benchleftflankseq = refobj.fetch(reference=chrom, start=int(start)-5, end=int(start))
            benchleftflankseq = benchleftflankseq.upper()
            benchrightflankseq = refobj.fetch(reference=chrom, start=int(end), end=int(end)+5)
            benchrightflankseq = benchrightflankseq.upper()
            runlength = int(end) - int(start)
            if runlength not in mononucdict:
                mononucdict[runlength] = {}
            namefields = name.split("_")
            repeatedbase = namefields[-1]
            refalleleseq = repeatedbase * runlength
            readsforthismononuc = 0
            # for each read, attempt to find read positions aligned to the bases immediately flanking
            # the mononucleotide. If one side and/or the other has no aligning base (due either to a 
            # deletion or to clipped sequence) then check for clipping beginning at that base, and
            # extend the mononuc-aligned read sequence with clipped sequence to include adjacent 
            # repeated bases plus 5 bp additional (this is NIST's so-called "slop" length).
            # If a flanking base *does* have a read base aligned to it, also extend that end to 
            # include repeated bases plut 5 bp slop.
            for readalign in align_obj.fetch(contig=chrom, start=int(start), stop=int(end)):
                if max_reads_per_mononuc > 0 and readsforthismononuc >= max_reads_per_mononuc:
                    break
                if args.downsample is not None and random.random() >= args.downsample:
                    continue
                # query_sequence includes soft clipped bases:
                queryseq = readalign.query_sequence
                queryseq = queryseq.upper()
                if readalign.is_secondary or queryseq is None:
                    continue
                readname = readalign.query_name
                # if the read doesn't align to the entire mononucleotide, skip it:
                if readalign.reference_start is None or readalign.reference_end is None:
                    logger.debug("Read " + readname + " has no reference_start or reference_end for " + chrom + ":" + str(start) + "-" + str(end))
                    continue
                if readalign.reference_start > int(start) or readalign.reference_end < int(end):
                    continue
                # if the read doesn't align to the flanking bases, it must be clipped:
                cigartuples = readalign.cigartuples
                if ((readalign.reference_start == int(start) and cigartuples[0][0] != 4) or
                    (readalign.reference_end == int(end) and cigartuples[len(cigartuples)-1][0] != 4)):
                    logger.debug("Read " + readname + " aligns just to end of " + name + " but is not soft clipped--skipping")
                    continue
                pairs = readalign.get_aligned_pairs()
                if len(pairs) == 0:
                    continue
                # find zero-based read pos of base aligned to ref base one before mononuc:
                if readalign.reference_start < int(start):
                    readstart = find_readpos_in_pairs(pairs, int(start)-1)
                elif readalign.query_alignment_start > 0: # must be clipped--set start to one base before reference start:
                    readstart = readalign.query_alignment_start - 1
                # find zero-based read pos of base aligned to ref base one after mononuc:
                if readalign.reference_end > int(end):
                    readend = find_readpos_in_pairs(pairs, int(end))
                elif readalign.query_alignment_end < len(queryseq):
                    readend = readalign.query_alignment_end + 1

                if readstart is not None and readend is not None:
                    alleleseq = queryseq[readstart+1:readend]
                    # extend the read sequence so 5' and 3' to include all of repeated bases
                    #if readstart >=0:
                        #logger.debug("Checking " + queryseq[readstart] + " vs " + repeatedbase)
                    while readstart >= 0 and queryseq[readstart] == repeatedbase:
                        #logger.debug("Extending " + alleleseq + " one base")
                        alleleseq = repeatedbase + alleleseq
                        readstart = readstart - 1
                    #if readend < len(queryseq):
                        #logger.debug("Checking " + queryseq[readend] + " vs " + repeatedbase)
                    while readend < len(queryseq) and queryseq[readend] == repeatedbase:
                        #logger.debug("Extending " + alleleseq + " one base")
                        alleleseq = alleleseq + repeatedbase
                        readend = readend + 1

                    # check that read has 5 flanking base pairs on either side of the HP:
                    if readstart < 5 or readend > len(queryseq)-5:
                        logger.debug("Skipping align of read " + readname + " which has fewer than 5 bp outside of HP repeat for " + chrom + ":" + str(start) + "-" + str(end))
                        continue
                    else:
                        leftflankseq = queryseq[readstart-4:readstart+1]
                        leftflankseq = leftflankseq.upper()
                        rightflankseq = queryseq[readend:readend+5]
                        rightflankseq = rightflankseq.upper()
                    matchtype = "UNKNOWN"
                    if p[repeatedbase].match(alleleseq):
                        numbases = len(alleleseq)
                        if numbases != runlength:
                            potentialhetname = chrom + "_" + str(int(start)+1) + "_" + refalleleseq + "_" + alleleseq
                            if potentialhetname in hetsitedict:
                                matchtype = "HET"
                            else:
                                matchtype = "LENGTHERROR"
                        else:
                            if leftflankseq == benchleftflankseq and rightflankseq == benchrightflankseq:
                                matchtype = "CORRECT"
                            else:
                                logger.debug(leftflankseq + "/" + rightflankseq + " doesnt match " + benchleftflankseq + "/" + benchrightflankseq)
                                matchtype = "FLANKERROR"
                    else: # discrepancies within the HP 
                        numbases = -1
                        matchtype = "COMPLEX"

                    msfh.write(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + readname + "\t" + repeatedbase + "\t" + str(runlength) + "\t" + str(numbases) + "\t" + matchtype + "\n")
                    if numbases not in mononucdict[runlength]:
                        mononucdict[runlength][numbases] = {}
                    if matchtype not in mononucdict[runlength][numbases]:
                        mononucdict[runlength][numbases][matchtype] = 1
                    else:
                        mononucdict[runlength][numbases][matchtype] = mononucdict[runlength][numbases][matchtype] + 1
                    readsforthismononuc = readsforthismononuc + 1
                else:
                    logger.debug(readname + " is unaligned at one endpoint of interval " + str(start) + " to " + str(end) + " and has fewer than 5 clipped bases!")
                    if readstart is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a start")
                    if readend is None:
                        logger.warning(readname + " " + str(readalign.reference_start) + "-" + str(readalign.reference_end) + " does not have a end")

            mononucline = mfh.readline()

    return mononucdict

def retrieve_stats_from_strstats_file(filename)->dict:
    strdict = {}

    with open(filename, "r") as sfh:
        statline = sfh.readline()
        while statline:
            statline = statline.rstrip()
            [chrom, start, end, readname, repeatedbases, runlength, numbases, matchtype] = statline.split("\t")
            start = int(start)
            end = int(end)
            runlength = int(runlength)
            if runlength not in strdict:
                strdict[runlength] = {}
            numbases = int(numbases)
            if numbases not in strdict[runlength]:
                strdict[runlength][numbases] = {}
            if matchtype not in strdict[runlength][numbases]:
                strdict[runlength][numbases][matchtype] = 1
            else:
                strdict[runlength][numbases][matchtype] = strdict[runlength][numbases][matchtype] + 1
            statline = sfh.readline()

    return strdict

# "pairs" is a structure created by pysam for an alignment, and contains a list
#   of tuples "consisting of the 0-based offset from the start of the read 
#   sequence followed by the 0-based reference position
#
# This routine searches the list of tuples for the one corresponding to the 
#   desired reference position ("pos"), and returns the read position
#   for that tuple, assuming it exists
#
# Why would a tuple have a ref position of "None"? These tuples are insertions
#   in the read sequence (and deletions from the reference have read pos of 
#   None). In the case of deletions from the reference, the "ifnone" argument
#   can tell the routine whether to report the read position that is one to 
#   the left of the deleted ref location ("lower") or one to the right ("higher")
#
# Note: this routine takes a *zero-based* position as its "pos" argument!
#
def find_readpos_in_pairs(pairs, pos, ifnone="lower"):
    # number of tuples:
    alignlength = len(pairs)
    ilow = 0
    ihigh = alignlength - 1

    hiref = pairs[ihigh][1]
    while hiref is None and ihigh > 0:
        ihigh = ihigh - 1
        hiref = pairs[ihigh][1]
    lowref = pairs[ilow][1]
    while lowref is None and ilow < len(pairs):
        ilow = ilow + 1
        lowref = pairs[ilow][1]

    # these shouldn't really happen
    if lowref is None or hiref is None or lowref > pos or hiref < pos or lowref > hiref:
        return None

    # need to be careful here to avoid an infinite loop:
    lastimid = None
    while ihigh - 1 > ilow and hiref >= pos and lowref <= pos:
        imid = int((ilow + ihigh)/2)
        midref = pairs[imid][1]
        while midref is None:
            if ifnone == "lower":
                imid = imid - 1
            else:
                imid = imid + 1
            midref = pairs[imid][1]
            if (ifnone == "lower" and imid <= ilow) or (ifnone=="higher" and imid >= ihigh):
                break
        if midref is None:
            return None
        if lastimid is not None and lastimid == imid:
            for i in range(ilow, ihigh):
                if pairs[i][1] is not None and pairs[i][1] == pos:
                    return pairs[i][0]
            return None
        
        if midref == pos:
            return pairs[imid][0] # might be None if nothing is aligned here
        elif midref > pos:
            ihigh = imid
            hiref = pairs[ihigh][1]
        elif midref < pos:
            ilow = imid
            lowref = pairs[ilow][1]
        if (pairs[ilow][1] is not None and pairs[ilow][1] > pos) or (pairs[ihigh][1] is not None and pairs[ihigh][1] < pos):
            return None
        lastimid = imid

    return None

def widen_str_on_readseq(readseq:str, readstart:int, readend:int, repbases:str)->list:
    alleleseq = readseq[readstart+1:readend]

    replength = len(repbases) 
    readlength = len(readseq) 
    perfectrepeat = True

    initialalleleseq = alleleseq
    initialreadstart = readstart
    initialreadend = readend

    # check whether read's current allele sequence is purely composed of the repeat:
    alleleindex = 0
    while alleleindex < len(alleleseq):
        allelebase = alleleseq[alleleindex]
        repbaseindex = alleleindex % replength
        repbase = repbases[repbaseindex]
        if repbase != allelebase:
            logger.debug("Allele sequence " + alleleseq + " is not a perfect repeat of " + repbases)
            perfectrepeat = False
            return [alleleseq, readstart, readend, perfectrepeat]
        alleleindex = alleleindex + 1

    # extend to left:
    readindex = readstart
    repbaseindex = (readindex - 1 - readstart) % replength
    while readindex > 0 and readseq[readindex] == repbases[repbaseindex]:
        readstart = readstart - 1
        alleleseq = readseq[readindex] + alleleseq
        readindex = readindex - 1
        repbaseindex = (readindex - 1 - initialreadstart) % replength

    # extend to right:
    readindex = readend
    repbaseindex = (readindex - readstart - 1) % replength
    while readindex < readlength and readseq[readindex] == repbases[repbaseindex]:
        readend = readend + 1
        alleleseq = alleleseq + readseq[readindex]
        readindex = readindex + 1
        repbaseindex = (readindex - readstart - 1) % replength

    if initialalleleseq != alleleseq:
        logger.debug("Extended " + initialalleleseq + " to " + alleleseq + " " + str(initialreadstart)  + "-" + str(initialreadend) + "/" + str(readstart) + "-" + str(readend) )

    return [alleleseq, readstart, readend, perfectrepeat]

def assess_read_align_errors(align_obj, refobj, readerrorfile:str, bedintervals, hetsitedict, args):
   
    if not args.errorfile and not args.rerun:
        refh = open(readerrorfile, "w")
    elif args.errorfile:
        readerrorfile = args.errorfile 

    stats = {}
    stats["totalalignedbases"] = 0
    stats["totalclippedbases"] = 0
    stats["totalerrorsinaligns"] = 0
    stats["singlebasecounts"] = {}
    stats["indellengthcounts"] = {}
    stats["alignedqualscorecounts"] = []
    # position counts are zero-based, i.e., one less than the ordinal position in the read:
    stats["positionsnvcounts"] = []
    stats["positionindelcounts"] = []
    stats["positiontotalcounts"] = []
    stats["snverrorqualscorecounts"] = []
    stats["indelerrorqualscorecounts"] = []
    alignsprocessed = 0

    # are we dealing with just a set of regions? or the whole genome?
    if bedintervals is None:
        bedintervals = [None]

    for benchinterval in bedintervals:
        if benchinterval is not None:
            regionstring = benchinterval.chrom + ":" + str(benchinterval.start + 1) + "-" + str(benchinterval.end)
        else:
            regionstring = None

        for align in align_obj.fetch(region=regionstring):

            if align.is_secondary or align.cigartuples is None or (args.downsample is not None and random.random() >= args.downsample):
                continue

            stats["totalalignedbases"] = stats["totalalignedbases"] + align.reference_length
            if benchinterval is not None: # shorten aligned length contribution if necessary
                if benchinterval.start > align.reference_start:
                    stats["totalalignedbases"] = stats["totalalignedbases"] + align.reference_start - benchinterval.start
                if benchinterval.end < align.reference_end:
                    stats["totalalignedbases"] = stats["totalalignedbases"] - align.reference_end + benchinterval.end

            # TODO: clipping calc should be amended to not count if alignment endpoints are outside desired benchinterval!
            cigartuples = align.cigartuples
            if cigartuples[0][0] in [4, 5]:
                if benchinterval is None or benchinterval.start < align.reference_start:
                    stats["totalclippedbases"] = stats["totalclippedbases"] + cigartuples[0][1]
            if cigartuples[-1][0] in [4, 5]:
                if benchinterval is None or benchinterval.end > align.reference_end:
                    stats["totalclippedbases"] = stats["totalclippedbases"] + cigartuples[-1][1]
            
            if not args.errorfile and not args.rerun:
                query, querystart, queryend, ref, refstart, refend, strand = alignparse.retrieve_align_data(align)
                if strand == "F":
                    queryleft = querystart
                    queryright = queryend
                else:
                    queryleft = queryend
                    queryright = querystart

                # make sure position arrays (zero-base indexed) are long enough to store position counts:
                while len(stats["positionindelcounts"]) < queryend:
                    stats["positionindelcounts"].append(0)
                while len(stats["positionsnvcounts"]) < queryend:
                    stats["positionsnvcounts"].append({})
                while len(stats["positiontotalcounts"]) < queryend:
                    stats["positiontotalcounts"].append(0)
                for position in range(querystart-1, queryend):
                    stats["positiontotalcounts"][position-1] = stats["positiontotalcounts"][position-1] + 1

                hetsites = {}
                hetsitealleles = {} # no need to track het site alleles in this context
                queryobj = None
    
                read_variants = alignparse.align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand, hetsites, hetsitealleles, stats["alignedqualscorecounts"], stats["snverrorqualscorecounts"], stats["indelerrorqualscorecounts"], True)
    
                for variant in read_variants:
                    namefields = variant.name.split("_")
                    # this is the read position (from beginning of read, regardless of strand):
                    pos = int(namefields[-4])
                    refallele = namefields[-3]
                    altallele = namefields[-2]
                    if benchinterval is not None and (variant.start < benchinterval.start or variant.end > benchinterval.end):
                        continue
                    if namefields[-1] == "F":
                        alignstrand = '+'
                    else:
                        alignstrand = '-'
                    varname = variant.chrom + "_" + str(int(variant.start) + 1) + "_" + refallele + "_" + altallele
                    varqvscore = "1000"
                    if variant.qvscore is not None:
                        varqvscore = str(variant.qvscore)

        
                    if varname in hetsitedict.keys():
                        errortype = 'HET'
                        errortypecolor = '255,0,0'
                    else:
                        errortype = 'ERROR'
                        errortypecolor = '0,0,255'
    
                    # record error in stats:
                    if errortype == 'ERROR':
                        stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                        if refallele != "*" and altallele != "*" and len(refallele) == 1 and len(altallele) == 1:
                            if alignstrand == "+":
                                readrefallele = refallele
                                readaltallele = altallele
                            else:
                                readrefallele = seqparse.revcomp(refallele)
                                readaltallele = seqparse.revcomp(altallele)

                            snvkey = readrefallele + "_" + readaltallele
                            if snvkey in stats["singlebasecounts"]:
                                stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                            else:
                                stats["singlebasecounts"][snvkey] = 1
                            # add SNV error to tally at this read position:
                            if snvkey in stats["positionsnvcounts"][pos - 1]:
                                stats["positionsnvcounts"][pos-1][snvkey] = stats["positionsnvcounts"][pos-1][snvkey] + 1
                            else:
                                stats["positionsnvcounts"][pos-1][snvkey] = 1
                        else:
                            trueref = refallele.replace('*', '')
                            truealt = altallele.replace('*', '')
                            lengthdiff = len(truealt) - len(trueref)
                            if lengthdiff in stats["indellengthcounts"]:
                                stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                            else:
                                stats["indellengthcounts"][lengthdiff] = 1
                            #stats["positionindelcounts"][pos-1] = stats["positionindelcounts][pos-1] + 1
    
                    refh.write(variant.chrom + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + varname + "\t" + varqvscore + "\t" + alignstrand + "\t" + str(variant.start) + "\t" + str(variant.end) + "\t" + errortypecolor + "\t" + errortype + "\t" + variant.name + "\n")
    
            alignsprocessed = alignsprocessed + 1
            if alignsprocessed == 100000*int(alignsprocessed/100000):
                logger.debug("Processed " + str(alignsprocessed) + " aligns")


    if args.errorfile or args.rerun:
        with open(readerrorfile, "r") as efh:
            errorline = efh.readline()
            while errorline:
                errorline = errorline.rstrip()
                [chrom, start, end, varname, score, strand, widestart, wideend, color, variant_type, queryvariantname] = errorline.split("\t")
                if variant_type == "HET":
                    errorline = efh.readline()
                    continue

                namefields = varname.split("_")
                refallele = namefields[-2]
                altallele = namefields[-1]
                
                # record error in stats:
                stats["totalerrorsinaligns"] = stats["totalerrorsinaligns"] + 1
                if refallele != "*" and altallele != "*" and len(refallele) == 1 and len(altallele) == 1:
                    snvkey = refallele + "_" + altallele
                    if snvkey in stats["singlebasecounts"]:
                        stats["singlebasecounts"][snvkey] = stats["singlebasecounts"][snvkey] + 1
                    else:
                        stats["singlebasecounts"][snvkey] = 1
                else:
                    trueref = refallele.replace('*', '')
                    truealt = altallele.replace('*', '')
                    lengthdiff = len(truealt) - len(trueref)
                    if lengthdiff in stats["indellengthcounts"]:
                        stats["indellengthcounts"][lengthdiff] = stats["indellengthcounts"][lengthdiff] + 1
                    else:
                        stats["indellengthcounts"][lengthdiff] = 1
                errorline = efh.readline()

    return stats
