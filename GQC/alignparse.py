import re
import os
import pysam
import pybedtools
import logging
import statistics
from collections import namedtuple
from pathlib import Path
from GQC import seqparse
from GQC import phasing
from GQC import bedtoolslib
from GQC import output

# create namedtuple for bed intervals:
varianttuple = namedtuple('varianttuple', ['chrom', 'start', 'end', 'name', 'vartype', 'excluded', 'qvscore']) 

logger = logging.getLogger(__name__)

# query start is always less than query end regardless of strand. query left always corresponds to ref start, 
# and so will be greater than query right when strand is reversed--all four are 1-based
def write_bedfiles(bamobj, pafaligns, refobj, queryobj, hetsites, testmatbed, testpatbed, truthbed, hetallelebed, excludedbedobj, args):

    refcoveredstring = ""
    querycoveredstring = ""
    refcoveredbed = None
    querycoveredbed = None
    if 'variantfile' in args:
        user_variantfile = args.variantfile
    else:
        user_variantfile = None
    variants = []
    hetsitealleles = {}
    alignedscorecounts = []
    snverrorscorecounts = []
    indelerrorscorecounts = []

    if ((testmatbed is not None and not os.path.exists(testmatbed)) or (testpatbed is not None and not os.path.exists(testpatbed)) or (truthbed is not None and not os.path.exists(truthbed))):
        if bamobj is not None:
            for align in bamobj.fetch():
                if align.is_secondary:
                    continue
                if align.reference_length >= args.minalignlength:
                    query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
                    if strand == "F":
                        queryleft = querystart
                        queryright = queryend
                    else:
                        queryleft = queryend
                        queryright = querystart
                    querynamestring = query + "." + str(queryleft) + "." + str(queryright)
                    refnamestring = ref + "." + str(refstart) + "." + str(refend) + "." + strand
                    querycoveredstring += query + "\t" + str(querystart - 1) + "\t" + str(queryend) + "\t" + refnamestring + "\n"
                    refcoveredstring += ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querynamestring + "\n"
                    if user_variantfile is None:
                        variants.extend(align_variants(align, queryobj, query, querystart, queryend, refobj, ref, refstart, refend, strand, hetsites, hetsitealleles, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts, True))
            # mark variants that are in excluded regions:
            logger.debug("Beginning to exclude variants in excluded regions")
            if excludedbedobj:
                variants = exclude_variants(variants, excludedbedobj)
                logger.debug("Finished excluding variants in excluded regions")
        else:
            for pafdict in pafaligns:
                # all start/endpoints are 1-based
                query = pafdict['query']
                querystart = pafdict['querystart']
                queryend = pafdict['queryend']
                if querystart < queryend:
                    queryleft = querystart
                    queryright = queryend
                else:
                    queryleft = queryend
                    queryright = querystart
                ref = pafdict['target']
                refstart = pafdict['targetstart']
                refend = pafdict['targetend']
                strand = pafdict['strand']

                querynamestring = query + "." + str(querystart) + "." + str(queryend)
                refnamestring = ref + "." + str(refstart) + "." + str(refend) + "." + strand
                querycoveredstring += query + "\t" + str(queryleft - 1) + "\t" + str(queryright) + "\t" + refnamestring + "\n"
                refcoveredstring += ref + "\t" + str(refstart - 1) + "\t" + str(refend) + "\t" + querynamestring + "\n"

        refcoveredbed = pybedtools.BedTool(refcoveredstring, from_string = True)
        querycoveredbed = pybedtools.BedTool(querycoveredstring, from_string = True)
    else:
        logger.info("Skipping alignment parsing because alignment coverage bedfiles already exist--delete them to reprocess")
        refcoveredbed = pybedtools.BedTool(truthbed)

    # if testpatbed is None, a single bed file should be written to the file passed as testmatbed:
    if not os.path.exists(testmatbed) or (testpatbed is not None and not os.path.exists(testpatbed)):
        phap1 = re.compile(r'.*MAT.*')
        phap2 = re.compile(r'.*PAT.*')
        with open(testmatbed, "w") as tmb:
            for testint in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
                if testpatbed is None or phap1.match(testint.name):
                    tmb.write(testint.chrom + "\t" + str(testint.start) + "\t" + str(testint.end) + "\t" + testint.name + "\n")

        if testpatbed is not None:
            with open(testpatbed, "w") as tpb:
                for testint in sorted(querycoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
                    if phap2.match(testint.name):
                        tpb.write(testint.chrom + "\t" + str(testint.start) + "\t" + str(testint.end) + "\t" + testint.name + "\n")

    if not os.path.exists(truthbed):
        with open(truthbed, "w") as rb:
            for truthint in sorted(refcoveredbed, key=lambda h: (h.chrom, h.start, h.stop)):
                rb.write(truthint.chrom + "\t" + str(truthint.start) + "\t" + str(truthint.end) + "\t" + truthint.name + "\n")

    if hetsites and hetallelebed is not None and not os.path.exists(hetallelebed):
        phasing.write_hetallele_bed(hetsitealleles, hetallelebed)

    if user_variantfile is not None:
        with open(user_variantfile, "r") as vh:
            variantline = vh.readline()
            while variantline:
                variantline = variantline.rstrip()
                [chrom, start, end, name] = variantline.split("\t")
                variants.append(varianttuple(chrom=chrom, start=int(start), end=int(end), name=name, qvscore=None))
                variantline = vh.readline()

    return [refcoveredbed, querycoveredbed, variants, hetsitealleles, alignedscorecounts, snverrorscorecounts, indelerrorscorecounts]

def retrieve_align_data(align)->list:
    if align.is_reverse:
        strand = 'R'
    else:
        strand = 'F'
    query = align.query_name

    # find number of hard clipped bases:
    cigartuples = align.cigartuples
    # will use actual one-based positions, so I don't go crazy, then report BED format zero-based half open later
    # pysam's align.query_alignment_start is the 0-based coordinate within the uncomplemented hard clipped query sequence, so here we add hardclipping from the 
    # appropriate end to get coordinates within the entire sequence (so low coordinates are towards the start of the original query sequence, not the left end
    # of the alignment)
    if strand == 'F':
        if cigartuples[0][0] == 5:
            hardclip = cigartuples[0][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_alignment_start + 1 # first unclipped base (was 0-based, now 1-based)
        queryend = hardclip + align.query_alignment_end # was 0-based index just past the last unclipped base now 1-based last unclipped base
    elif strand == 'R':
        if cigartuples[-1][0] == 5:
            hardclip = cigartuples[-1][1]
        else:
            hardclip = 0
        querystart = hardclip + align.query_length - align.query_alignment_end + 1 # was 0-based index just past end of reversed query, now 1-based 5' end of query
        queryend = hardclip + align.query_length - align.query_alignment_start # was 0-based index of first unclipped base, now 1-based 3' unclipped end of query

    ref = align.reference_name

    refstart = align.reference_start + 1 # was zero-based, now one-based
    refend = align.reference_end #  was zero-based, adding one for bed format "half open", now one-based

    aligndata = [query, querystart, queryend, ref, refstart, refend, strand]

    return aligndata

# query start and query end are the lower and higher endpoints of the query seq in query coordinates (1-based)
# regardless of orientation of the alignment
def align_variants(align, queryobj, query:str, querystart:int, queryend:int, refobj, ref:str, refstart:int, refend:int, strand:str, chromhetsites={}, hetsitealleles={}, alignedscorecounts=[], snverrorscorecounts=[], indelerrorscorecounts=[], widen=True)->list:

    # coordinates are all one-based, with start at beginning of *original* sequence (not left end of the alignment)
    variantlist = []
    coveredregionlist = []
    homozygousregionlist = []

    # make an array of query positions for each ref position:
    query_positions = []

    if queryobj is None:
        queryseq = align.query_alignment_sequence
    else:
        queryseq = queryobj.fetch(reference=query, start=querystart-1, end=queryend).upper()

    # qual scores might be None, but if they're not, add to the counts in alignedscorecounts and errorscorecounts:
    alignedqualscores = align.query_alignment_qualities

    if alignedqualscores is not None:
        if len(alignedscorecounts)==0:
            for i in range(100):
                alignedscorecounts.append(0)
                snverrorscorecounts.append(0)
                indelerrorscorecounts.append(0)
        for qualscorechar in alignedqualscores:
            qualscore = int(qualscorechar)
            alignedscorecounts[qualscore] = alignedscorecounts[qualscore] + 1

    refseq = refobj.fetch(reference=ref, start=refstart-1, end=refend).upper()

    strandsign = 1
    bedstrand = '+'
    if strand == 'R':
        if queryobj is not None:
            queryseq = seqparse.revcomp(queryseq)
        strandsign = -1
        bedstrand = '-'

    alignops = align.cigartuples

    # first position will begin at left-most ref/query base of the alignment (regardless of strand)
    refcurrentoffset = 0
    querycurrentoffset = 0

    alignopindex = 0
    refalignlength = refend - refstart + 1
    queryalignlength = queryend - querystart + 1
    matchns = re.compile(".*[nN].*")
    numops = len(alignops)

    if query=="h2tg000004l":
        manualdebug = True
    else:
        manualdebug = False
    if manualdebug:
        print("Query has length " + str(queryalignlength) + " ref has length " + str(refalignlength))

    logger.debug("Traversing alignment with " + str(numops) + " cigar ops aligning " + ref + " to " + query)
    while refcurrentoffset <= refend-refstart and alignopindex < len(alignops): # traverse the alignment operator by operator
        alignop = alignops[alignopindex]
        op = alignop[0]
        oplength = alignop[1]

        if op in [0, 7, 8]: # MX= find SNV and MNVs
            if manualdebug:
                print("Pos " + ref + ":" + str(refcurrentoffset + refstart - 1) + " op M " + str(op) + " has length " + str(oplength))
            for blockoffset in range(oplength):
                refpos = refcurrentoffset + blockoffset # this is distance from left-most base of the alignment
                querypos = querycurrentoffset + blockoffset # this is distance from left-most base of the alignment
                if strand == 'F':
                    querycoordinate = querypos + querystart
                else:
                    querycoordinate = queryend - querypos
                if refpos >= len(refseq):
                    logger.debug("Ref position " + str(refpos) + " in " + ref + " is past end of seq")
                    continue
                if querypos >= len(queryseq):
                    logger.debug("Query position " + str(querypos) + " in " + query + " is past end of seq")
                    continue
                if refseq[refpos] != queryseq[querypos] and refseq[refpos] != "N" and queryseq[querypos] != "N":
                    variantname=query+"_"+str(querycoordinate)+"_"+refseq[refpos]+"_"+queryseq[querypos]+"_"+strand # query's 1-based position, ref base, query base (comp if rev strand), strand
                    #additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart-1) + "\t" + str(refpos+refstart) + "\t0,0,0\t" + alignstring
                    if alignedqualscores is not None:
                        snverrorscore = int(alignedqualscores[querypos])
                        snverrorscorecounts[snverrorscore] = snverrorscorecounts[snverrorscore] + 1
                    else:
                        snverrorscore = None
                    variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-1, end=refpos+refstart, name=variantname, vartype='SNV', excluded=False, qvscore=snverrorscore))
                elif queryseq[querypos] == 'N':
                    if alignedqualscores is not None:
                        snverrorscore = int(alignedqualscores[querypos])
                        snverrorscorecounts[snverrorscore] = snverrorscorecounts[snverrorscore] + 1

                query_positions.append(querypos)

        if op in [2, 3]: # deletions
            if manualdebug:
                print("Pos " + ref + ":" + str(refcurrentoffset + refstart - 1) + " op D " + str(op) + " has length " + str(oplength))
            refpos = refcurrentoffset-1;
            for deloffset in range(oplength):
                query_positions.append(querycurrentoffset)
            refallele = refseq[refcurrentoffset:refcurrentoffset+oplength] # one-based refstart+refcurrentoffset to refstart+refcurrentoffset+oplength-1
            queryallele = "*" # one-based between querystart+querycurrentoffset-1 and querystart+querycurrentoffset if forward strand, queryend-querycurrentoffset+1 and queryend-querycurrentoffset if rev
            extendleft = 0
            extendright = 0
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extendright" and *lower* the lefthand coordinate of ref entries and forward strand queries by "extendleft"
                while querycurrentoffset + extendright < queryalignlength and refcurrentoffset + extendright < refalignlength and refseq[refcurrentoffset + extendright] == queryseq[querycurrentoffset + extendright]:
                    refallele = refallele + queryseq[querycurrentoffset + extendright]
                    if queryallele == "*":
                        queryallele = queryseq[querycurrentoffset + extendright]
                    else:
                        queryallele = queryallele + queryseq[querycurrentoffset + extendright]
                    extendright = extendright + 1
                while querycurrentoffset - 1 - extendleft >= 0 and refcurrentoffset - 1 + oplength - extendleft >= 0 and refseq[refcurrentoffset - 1 + oplength - extendleft] == queryseq[querycurrentoffset - 1 - extendleft]:
                    refallele = queryseq[querycurrentoffset - extendleft - 1] + refallele
                    if queryallele == "*":
                        queryallele = queryseq[querycurrentoffset - extendleft - 1]
                    else:
                        queryallele = queryseq[querycurrentoffset - extendleft - 1] + queryallele
                    extendleft = extendleft + 1
            
            # querycoordinate values still need to be adjusted due to widening?
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset - extendleft
            else:
                querycoordinate = queryend - querycurrentoffset - extendright

            # if there are quality scores, find the median quality across this deletion for our error tally (and to report in the output)
            indelerrorscore = None
            if alignedqualscores is not None:
                if extendleft == 0 and extendright == 0:
                    queryquals = alignedqualscores[querycurrentoffset-1:querycurrentoffset]
                    qualscores = [int(x) for x in queryquals]
                else:
                    queryquals = alignedqualscores[querycurrentoffset-extendleft:querycurrentoffset+extendright]
                    qualscores = [int(x) for x in queryquals]
                if len(qualscores)==0:
                    logger.debug("Variant with pos " + ref + ":" + str(refpos+refstart-extendleft) + "-" + str(refpos+refstart+oplength+extendright) + " name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has query surrounding seq " + querysurroundingseq + " and length of qualscores is zero!")
                else:
                    numquals = len(qualscores)
                    # if even number of qual scores, drop the top one so the lower of the two medians is chosen (rather than an average, which may not be represented in the total qv score counts)
                    if numquals==2*int(numquals/2):
                        qualscores.pop()
                    medqual = int(statistics.median(qualscores))
                    indelerrorscorecounts[medqual] = indelerrorscorecounts[medqual] + 1
                    indelerrorscore = medqual

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            # (looks as if I might be looking only at left and right bases within the query alleles here--need to check)
            [queryleftbase, queryrightbase] = ["", ""]
            if querycurrentoffset > 0:
                queryleftbase = queryseq[querycurrentoffset-1]
            if querycurrentoffset+extendright < queryalignlength:
                queryrightbase = queryseq[querycurrentoffset+extendright]
            querysurroundingseq = queryleftbase + queryrightbase
        
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(querysurroundingseq)):
                variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand # positions of insertions are positions to the left of first inserted base

                variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-extendleft, end=refpos+refstart+oplength+extendright, name=variantname, vartype='INDEL', excluded=False, qvscore=indelerrorscore ))
            else:
                variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand # positions of insertions are positions to the left of first inserted base
                logger.debug("Variant with name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has query surrounding seq " + querysurroundingseq + " was excluded")

        if op == 1: # insertion
            if manualdebug:
                print("Pos " + ref + ":" + str(refcurrentoffset + refstart - 1) + " op I " + str(op) + " has length " + str(oplength))
            refpos = refcurrentoffset-1;
            refallele = "*"
            queryallele = queryseq[querycurrentoffset:querycurrentoffset+oplength] # one-based querystart+querycurrentoffset to querystart+querycurrentoffset+oplength-1 if forward strand, queryend-querycurrentoffset to queryend-querycurrentoffset-oplength+1 if reverse
            if alignedqualscores is not None:
                queryquals = alignedqualscores[querycurrentoffset:querycurrentoffset+oplength]
                qualscores = [int(x) for x in queryquals]
            #query_positions.append(querycurrentoffset)
            extendright = 0
            extendleft = 0
            if widen is True: # n.b. - this will *lower* the righthand coordinate of reverse strand queries by "extendright"
                while querycurrentoffset + extendright < queryalignlength and refcurrentoffset + extendright < refalignlength and refseq[refcurrentoffset + extendright] == queryseq[querycurrentoffset + extendright]:
                    queryallele = queryallele + refseq[refcurrentoffset + extendright]
                    if alignedqualscores is not None:
                        qualscores.append(int(alignedqualscores[querycurrentoffset + extendright]))
                    if refallele == "*":
                        refallele = refseq[refcurrentoffset + extendright]
                    else:
                        refallele = refallele + refseq[refcurrentoffset + extendright]
                    extendright = extendright + 1
                while querycurrentoffset - 1 + oplength - extendleft >= 0 and refcurrentoffset - 1 - extendleft >= 0 and refseq[refcurrentoffset - 1 - extendleft] == queryseq[querycurrentoffset - 1 + oplength - extendleft]:
                    queryallele = refseq[refcurrentoffset - extendleft - 1] + queryallele
                    if alignedqualscores is not None:
                        qualscores.insert(0, int(alignedqualscores[querycurrentoffset - 1 + oplength - extendleft]))
                    if refallele == "*":
                        refallele = refseq[refcurrentoffset - extendleft - 1]
                    else:
                        refallele = refseq[refcurrentoffset - extendleft - 1] + refallele
                    extendleft = extendleft + 1
            if strand == 'F':
                querycoordinate = querystart + querycurrentoffset - extendleft
                querycoordend = querycoordinate + oplength - 1 # this is potentially off by one and could be a bug (see its use below)
            else:
                querycoordinate = queryend - querycurrentoffset - extendleft
                querycoordend = querycoordinate - oplength - 1 # this is potentially off by one and could be a bug (see its use below)

            # if there are quality scores, find the median quality across this insertion for our error tally (and to report in the output)
            indelerrorscore = None
            if alignedqualscores is not None:
                numquals = len(qualscores)
                # if even number of qual scores, drop the top one so the lower of the two medians is chosen (rather than an average, which may not be represented in the total qv score counts)
                if numquals==2*int(numquals/2):
                    qualscores.pop()
                medqual = int(statistics.median(qualscores))
                indelerrorscorecounts[medqual] = indelerrorscorecounts[medqual] + 1
                indelerrorscore = medqual

            # check neighboring bases for Ns, which can create misleading/possibly wrong variants
            [refleftbase, refrightbase] = ["", ""]
            if refcurrentoffset-extendleft > 0:
                refleftbase = refseq[refcurrentoffset-1-extendleft]
            else:
                refleftbase = ''
            if refcurrentoffset+extendright+1 < refalignlength:
                refrightbase = refseq[refcurrentoffset+extendright+1]
            else:
                refrightbase = ''
            refsurroundingseq = refleftbase + refrightbase
            
            if not (matchns.match(queryallele) or matchns.match(refallele) or matchns.match(refsurroundingseq)):
                variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand
                #additionalfields = "0\t" + bedstrand + "\t" + str(refpos+refstart) + "\t" + str(refpos+refstart+extendright) + "\t0,0,0\t" + alignstring
                logger.debug("Variant with pos " + ref + ":" + str(refpos+refstart-extendleft) + "-" + str(refpos+refstart+extendright) + " name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has ref surrounding seq " + refsurroundingseq)
                variantlist.append(varianttuple(chrom=ref, start=refpos+refstart-extendleft, end=refpos+refstart+extendright, name=variantname, vartype='INDEL', excluded=False, qvscore=indelerrorscore ))
            else:
                variantname=query+"_"+str(querycoordinate)+"_"+refallele+"_"+queryallele+"_"+strand
                logger.debug("Variant with name " + variantname + " and queryallele " + queryallele + " and refallele " + refallele + " has ref surrounding seq " + refsurroundingseq + " was excluded")

        # advance current positions: cases where reference coord advances (MDN=X):
        if op in [0, 2, 3, 7, 8]:
            refcurrentoffset = refcurrentoffset + oplength
        # cases where query coord advances (MI=X)
        if op in [0, 1, 7, 8]:
            querycurrentoffset = querycurrentoffset + oplength

        alignopindex = alignopindex + 1

    numquerypositions = len(query_positions)
    refseqlength = len(refseq)

    # use query positions to assess het alleles/covered regions:
    if chromhetsites and ref in chromhetsites:
        desiredhets = chromhetsites[ref]
        for het in desiredhets:
            hetname = het.name
            namefields = hetname.split("_")
            hetpos = int(namefields[-4])
            refallele = namefields[-3]
            altallele = namefields[-2]

            hetstart = het.start
            hetend = het.end
            if hetstart < refstart:
                continue
            if hetend >= refend:
                break
            querystartoffset = query_positions[hetstart-refstart]
            if hetend - refstart + 2 >= len(query_positions):
                queryendoffset = len(queryseq)
            else:
                queryendoffset = query_positions[hetend-refstart+2]
            queryallele = '*'
            if querystartoffset + 1 < queryendoffset - 1:
                queryallele = queryseq[querystartoffset+1:queryendoffset-1]
            alleletype = 'NEITHER'
            if queryallele == refallele:
                alleletype = 'SAME'
            elif queryallele == altallele:
                alleletype = 'ALT'

            if strand == "F":
                querystartcoord = querystart + querystartoffset
                queryendcoord = querystart + queryendoffset - 2
            else:
                querystartcoord = queryend - queryendoffset
                queryendcoord = queryend - querystartoffset

            # adjust end coordinate for cases where there's a deletion in the test assembly at the position of the het:
            if queryendcoord < querystartcoord:
                queryendcoord = querystartcoord

            hetsitealleles[hetname] = {'name':hetname, 'ref':ref, 'refstart':hetstart, 'refend':hetend, 'allele':queryallele, 'query':query, 'start':querystartcoord, 'end':queryendcoord, 'chrom':query}

    return variantlist

def split_aligns_and_sort(splitbamname, bamobj, minindelsize=10000):
    with pysam.AlignmentFile(splitbamname, "wb", header=bamobj.header) as sbfh:
        for align in bamobj.fetch():
            if align.is_secondary:
                continue
            splitaligns = split_align_on_indels(align, minindelsize)
            numaligns = len(splitaligns)
            for splitalign in splitaligns:
                sbfh.write(splitalign)

def split_align_on_indels(align, minindelsize=10000)->list:

    subaligninfo = [] # array of data regarding the sub-alignments

    # general stats about this alignment:
    if align.is_reverse:
        strand = 'R'
    else:
        strand = 'F'
    query = align.query_name
    ref = align.reference_name
    refstart = align.reference_start + 1 # was zero-based, now one-based
    refend = align.reference_end #  was zero-based, adding one for bed format "half open", now one-based
    [left_softclip, right_softclip] = left_right_soft_clip(align)
    [left_hardclip, right_hardclip] = left_right_hard_clip(align)
    originalflag = align.flag
    alignops = align.cigartuples

    logger.debug("Alignment " + align.query_name + " to " + ref + " with soft clip " + str(left_softclip) + "/" + str(right_softclip) + " and hard clip " + str(left_hardclip) + "/" + str(right_hardclip) + " and flag " + str(originalflag))

    # will we hard or soft clip the longest alignment?
    hardcliplongest = False
    if align.is_supplementary or (left_hardclip > 0) or (right_hardclip > 0):
        hardcliplongest = True

    # first position will begin at left-most ref/query base of the alignment (regardless of strand)
    refcurrentoffset = 0
    querycurrentoffset = 0
    reflastalignend = 0 # also zero-based index for the next align segment--this is an offset from refstart
    querylastalignend = 0 # also zero-based index for the next align segment--this is an offset from refstart
    lastalignopindex = 0 # also zero-based index for the next align segment--this is an offset from the first cigar operator
    if left_softclip > 0 or left_hardclip > 0:
        lastalignopindex = 1
    currentalignopindex = 0

    segnumber = 1
    while refcurrentoffset <= refend-refstart and currentalignopindex < len(alignops): # traverse the alignment operator by operator
        alignop = alignops[currentalignopindex]
        op = alignop[0]
        oplength = alignop[1]

        if op in [2, 3]: # deletions
            # if deletion is large enough, append the portion of the alignment to the left as a sub-align, but without the deletion
            if oplength >= minindelsize:
                logger.debug("Align segment " + ref + ":" + str(refstart + reflastalignend) + "-" + str(refstart + refcurrentoffset) + " segment " + str(segnumber))
                logger.debug("Deletion size " + str(op) + " is larger than " + str(minindelsize) + ": splitting")
                segnumber = segnumber + 1
                subalignlength = refcurrentoffset - reflastalignend
                subaligninfo.append({'alignedquerystart':querylastalignend, 'alignedqueryend':querycurrentoffset, 'alignedrefstart':refstart+reflastalignend, 'alignedrefend':refstart+refcurrentoffset, 'cigarops':alignops[lastalignopindex:currentalignopindex], 'subalignlength':subalignlength, 'segnum':segnumber })
                # also need to advance the start coordinate for the next align to point beyond deletion:
                lastalignopindex = currentalignopindex + 1
                reflastalignend = refcurrentoffset + oplength
                querylastalignend = querycurrentoffset
                logger.debug("Setting querylastalignend to " + str(querylastalignend) + " after deletion of length " + str(oplength))

        if op == 1: # insertion
            # if insertion is large enough, append the portion of the alignment to the left as a split align, but without the insertion
            if oplength >= minindelsize:
                logger.debug("Align segment " + ref + ":" + str(refstart + reflastalignend) + "-" + str(refstart + refcurrentoffset) + " to " + str(querylastalignend) + "-" + str(querycurrentoffset) + " segment " + str(segnumber))
                logger.debug("Insertion size " + str(op) + " is larger than " + str(minindelsize) + ": splitting")
                segnumber = segnumber + 1
                subalignlength = refcurrentoffset - reflastalignend
                subaligninfo.append({'alignedquerystart':querylastalignend, 'alignedqueryend':querycurrentoffset, 'alignedrefstart':refstart+reflastalignend, 'alignedrefend':refstart+refcurrentoffset, 'cigarops':alignops[lastalignopindex:currentalignopindex], 'subalignlength':subalignlength, 'segnum':segnumber })
                # also need to advance the start coordinate for the next align to point beyond insertion:
                lastalignopindex = currentalignopindex + 1
                reflastalignend = refcurrentoffset
                querylastalignend = querycurrentoffset + oplength
                logger.debug("Setting querylastalignend to " + str(querylastalignend) + " after insertion of length " + str(oplength))
        currentalignopindex = currentalignopindex + 1

        # advance current positions: cases where reference coord advances (MDN=X):
        if op in [0, 2, 3, 7, 8]:
            refcurrentoffset = refcurrentoffset + oplength
        # cases where query coord advances (MI=X)
        if op in [0, 1, 7, 8]:
            querycurrentoffset = querycurrentoffset + oplength

    logger.debug("Align segment " + ref + ":" + str(refstart + reflastalignend) + "-" + str(refstart + refcurrentoffset) + " to " + str(querylastalignend) + "-" + str(querycurrentoffset) + " last segment " + str(segnumber))
    subalignlength = refcurrentoffset - reflastalignend
    subaligninfo.append({'alignedquerystart':querylastalignend, 'alignedqueryend':querycurrentoffset, 'alignedrefstart':refstart+reflastalignend, 'alignedrefend':refstart+refcurrentoffset, 'cigarops':alignops[lastalignopindex:currentalignopindex], 'subalignlength':subalignlength, 'segnum':segnumber })

    subaligns = create_subalignobjects(align, subaligninfo, hardcliplongest)

    return subaligns

def fix_adjacent_insertions_deletions(cigarops)->list:
    newcigarops = []
    skipnext = False
    # wherever same-size insertions and deletions are adjacent, combine them:
    for i in range(len(cigarops)):
        if skipnext:
            skipnext = False
            continue
        cigartuple = cigarops[i]
        opnum = cigartuple[0]
        oplength = cigartuple[1]
        # if an insertion and the following tuple is a deletion:
        if i < len(cigarops) - 1 and opnum == 1 and cigarops[i+1][0] == 2:
            nextoplength = cigarops[i+1][1]
            skipnext = True
            logger.debug("Combining cigar " + str(opnum) + " length " + str(oplength) + " with cigar " + str(cigarops[i+1][0]) + " length " + str(nextoplength))
            print("Combining cigar " + str(opnum) + " length " + str(oplength) + " with cigar " + str(cigarops[i+1][0]) + " length " + str(nextoplength))
            if nextoplength == oplength:
                newcigarops.append([0, nextoplength])
            elif nextoplength > oplength:
                newcigarops.append([0, oplength])
                newcigarops.append([2, nextoplength - oplength])
            else:
                newcigarops.append([0, nextoplength])
                newcigarops.append([1, oplength - nextoplength])
        # if a deletion and the following tuple is an insertion:
        elif i < len(cigarops) - 1 and opnum == 2 and cigarops[i+1][0] == 1:
            nextoplength = cigarops[i+1][1]
            skipnext = True;
            logger.debug("Combining cigar " + str(opnum) + " length " + str(oplength) + " with cigar " + str(cigarops[i+1][0]) + " length " + str(nextoplength))
            print("Combining cigar " + str(opnum) + " length " + str(oplength) + " with cigar " + str(cigarops[i+1][0]) + " length " + str(nextoplength))
            if nextoplength == oplength:
                newcigarops.append([0, nextoplength])
            elif nextoplength > oplength:
                newcigarops.append([0, oplength])
                newcigarops.append([1, nextoplength - oplength])
            else:
                newcigarops.append([0, nextoplength])
                newcigarops.append([2, oplength - nextoplength])
        # if not adjacent insertion/deletions, simply add this operation
        else:
            newcigarops.append([opnum, oplength])

    return newcigarops

def count_consumed_query(cigarops)->int:

    queryconsumed = 0
    for cigartuple in cigarops:
        opnum = cigartuple[0]
        oplength = cigartuple[1]
        if opnum in [0, 1, 4, 7, 8]:
            queryconsumed = queryconsumed + oplength
    return queryconsumed

def create_subalignobjects(align, subaligninfo:list, hardcliplongest:bool)->list:

    querysamseq = align.query_sequence
    samseqlength = len(querysamseq)
    [full_left_hardclip, full_right_hardclip] = left_right_hard_clip(align)
    subalignobjs = []
    longest = True
    cigarops = align.cigartuples
    for aligninfo in sorted(subaligninfo, key=lambda d: d['subalignlength'], reverse=True):
        newalign = pysam.AlignedSegment()
        newalign.query_name = align.query_name
        newalign.reference_id = align.reference_id
        newalign.mapping_quality = align.mapping_quality
        newalign.reference_start = aligninfo['alignedrefstart'] - 1
        cigartuples = fix_adjacent_insertions_deletions(aligninfo['cigarops'])
        oldqueryconsumed = count_consumed_query(aligninfo['cigarops'])
        newqueryconsumed = count_consumed_query(cigartuples)
        if oldqueryconsumed != newqueryconsumed:
            logger.info("Cigar op length change for " + align.query_name)
            logger.info("Old length: " + str(oldqueryconsumed) + "\nNew length: " + str(newqueryconsumed))

        # soft clip the longest unless hardcliplongest was specified
        if longest and not hardcliplongest:
            logger.debug("Will softclip align of length " + str(aligninfo['subalignlength']) + " with start " + str(aligninfo['alignedquerystart']))
            newalign.query_sequence = querysamseq
            leftsoftclip = align.query_alignment_start + aligninfo['alignedquerystart']
            rightsoftclip = samseqlength - aligninfo['alignedqueryend'] - align.query_alignment_start
            logger.debug("Will add " + str(leftsoftclip) + " left softclip and right softclip " + str(rightsoftclip) + " to align of length " +str(aligninfo['subalignlength']) + " with sam length " + str(samseqlength))
            if leftsoftclip > 0:
                cigarlist = [(4, leftsoftclip)] + cigartuples
            else:
                cigarlist = cigartuples
            if rightsoftclip > 0:
                cigarlist.append([4, rightsoftclip])
            try:
                newalign.cigartuples = cigarlist
            except OverflowError:
                logger.critical("Overflow Error while retrieving supplementary CIGAR string. Read name: {0}, Position: {1}, CIGAR: {2}".format(align.query_name, str(aligninfo['alignedquerystart']), str(aligninfo['alignedrefstart'])))
                continue
            #for cigar in cigarlist:
                #logger.debug("CIGAR\t" + str(cigar[0]) + ":" + str(cigar[1]))
            queryconsumedbases = count_consumed_query(cigarlist)
            logger.debug("Query seq length is " + str(samseqlength) + " and cigar consumes " + str(queryconsumedbases))
        else:
            logger.debug("Will hardclip align of length " + str(aligninfo['subalignlength']))
            lefthardclip = align.query_alignment_start + aligninfo['alignedquerystart'] + full_left_hardclip
            righthardclip = samseqlength - aligninfo['alignedqueryend'] - align.query_alignment_start + full_right_hardclip
            logger.debug("Will add " + str(lefthardclip) + " left hardclip and right hardclip " + str(righthardclip) + " to align of length " +str(aligninfo['subalignlength']) + " with sam length " + str(samseqlength))
            if lefthardclip > 0:
                cigarlist = [(5, lefthardclip)] + cigartuples
            else:
                cigarlist = cigartuples
            if righthardclip > 0:
                cigarlist.append([5, righthardclip])
            try:
                newalign.cigartuples = cigarlist
            except OverflowError:
                logger.critical("Overflow Error while retrieving supplementary CIGAR string. Read name: {0}, Position: {1}, CIGAR: {2}".format(align.query_name, str(aligninfo['alignedquerystart']), str(aligninfo['alignedrefstart'])))
                continue
            query_seqstart = aligninfo['alignedquerystart'] + align.query_alignment_start
            query_seqend = aligninfo['alignedqueryend'] + align.query_alignment_start
            newalign.query_sequence = querysamseq[query_seqstart:query_seqend]
            queryseqlength = query_seqend - query_seqstart
            #logger.debug("Query sequence has length " + str(lefthardclip) + " hardclip and right hardclip " + str(righthardclip) + " to align of length " +str(aligninfo['subalignlength']))
            queryconsumedbases = count_consumed_query(cigarlist)
            logger.debug("Query seq length is " + str(queryseqlength) + " and cigar consumes " + str(queryconsumedbases))

        if not longest:
            newalign.flag = align.flag | 2048
            logger.debug("Marking align of length " + str(aligninfo['subalignlength']) + " as supplementary")
        else:
            newalign.flag = align.flag
        longest = False
        subalignobjs.append(newalign)

    return subalignobjs

def exclude_variants(variants:list, excludedregionsobj:pybedtools.BedTool)->list:
    # create bedintervals for variants:
    logger.debug("Excluding lots of variants in excluded regions")
    numvariants = len(variants)
    logger.debug("There are " + str(numvariants) + " variants")
    variantbedstringlist = []
    for variant in variants:
        chrom = variant.chrom
        start = variant.start
        end = variant.end
        name = variant.name
        variantbedstringlist.append(chrom + "\t" + str(start) + "\t" + str(end) + "\t" + name + "\n")

    variantbedstring = ''.join(variantbedstringlist)
    variantbedobj = pybedtools.BedTool(variantbedstring, from_string=True)
    logger.debug("Created BedTool object")
    excludedvariants = bedtoolslib.intersectintervals(variantbedobj, excludedregionsobj, wa=True)
    excludeddict = {}
    for excludedvariant in excludedvariants:
        excludeddict[excludedvariant.name] = True

    newvariants = []
    numexcluded = 0
    numretained = 0
    for variant in variants:
        if variant.name in excludeddict.keys():
            newvariants.append(varianttuple(chrom=variant.chrom, start=variant.start, end=variant.end, name=variant.name, vartype=variant.vartype, excluded = True, qvscore=variant.qvscore))
            numexcluded = numexcluded + 1
        else:
            newvariants.append(variant)
            numretained = numretained + 1

    logger.debug("Excluded " + str(numexcluded) + " variants and kept " + str(numretained) + " variants")

    return newvariants

def read_paf_aligns(paffile:str, mintargetlength=0)->list:

    alignlist = []
    
    paf = Path(paffile)
    if not paf.is_file():
        logger.critical("PAF file " + paffile + " must exist and be readable")
        exit(1)

    with open(paffile, "r") as pfh:
        alignline = pfh.readline()
        while alignline:
            alignline = alignline.rstrip()
            fields = alignline.split("\t")
            if len(fields) >= 12:
                [query, querylength, querystartzb, queryend, strand, target, targetlength, targetstartzb, targetend, resmatches, blocklength, mapqual] = fields[0:12]
            else:
                logger.critical("Input paf-file has fewer than 12 tab-delimited columns. Unable to process.")
                exit(1)

            #if len(fields) > 12 and len(list(filter(lambda x:'tp:A:S' in x, fields[12:]))) > 0:
                #continue

            querystart = int(querystartzb) + 1
            queryend = int(queryend)
            queryalignlength = queryend - int(querystartzb)
            targetstart = int(targetstartzb) + 1
            targetend = int(targetend)
            targetalignlength = targetend - int(targetstartzb)

            if targetalignlength >= mintargetlength:
                if strand == "+":
                    alignlist.append({'query':query, 'querylength':int(querylength), 'querystart':querystart, 'queryend':queryend, 'queryalignlength':queryalignlength, 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':targetstart, 'targetend':targetend, 'targetalignlength':targetalignlength})
                else:
                    alignlist.append({'query':query, 'querylength':int(querylength), 'queryend':querystart, 'querystart':queryend, 'queryalignlength':queryalignlength, 'strand':strand, 'target':target, 'targetlength':int(targetlength), 'targetstart':targetstart, 'targetend':targetend, 'targetalignlength':targetalignlength})
            alignline = pfh.readline()

    return alignlist

def read_bam_aligns(bamobj, mintargetlength=0)->list:

    alignlist = []

    refentries = bamobj.references
    reflengths = bamobj.lengths
    reflengthdict = {}
    for tid in range(len(refentries)):
        reflengthdict[refentries[tid]] = reflengths[tid]

    for align in bamobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= mintargetlength:
            query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
            querylength = align.query_length
            targetlength = reflengthdict[ref]
            if strand == "F":
                queryleft = querystart
                queryright = queryend
                alignlist.append({'query':query, 'querylength':querylength, 'querystart':querystart, 'queryend':queryend, 'queryalignlength':queryend-querystart+1, 'strand':'+', 'target':ref, 'targetlength':targetlength, 'targetstart':refstart, 'targetend':refend, 'targetalignlength':refend-refstart+1})
            else:
                queryleft = queryend
                queryright = querystart
                alignlist.append({'query':query, 'querylength':querylength, 'queryend':querystart, 'querystart':queryend, 'queryalignlength':queryend-querystart+1, 'strand':'-', 'target':ref, 'targetlength':targetlength, 'targetstart':refstart, 'targetend':refend, 'targetalignlength':refend-refstart+1})

    return alignlist

# this routine assumes query start > query end for reverse strand alignments
def assess_overall_structure(aligndata:list, refobj, queryobj, outputfiles, bedobjects, benchmark_stats, args):

    # create a "clustercoverage" dictionary that will contain each chromosome's number of clusters and bases covered
    benchmark_stats["clustercoverage"] = {}
    benchmark_stats["alignclusters"] = {}
    benchmark_stats["numnonexcludedbases"] = {}

    # maximum separating distance along target to include in one cluster of alignments
    maxdistance = args.maxclusterdistance

    # create a directory for alignment line plots for each chromosome:
    output.create_output_directory(outputfiles["alignplotdir"])
    alignplotprefix = outputfiles["alignplotprefix"]

    # calculate number of non-excluded bases for each ref entry:
    # since allexcludedregions bed object is merged,can subtract out excluded bases interval by interval
    for ref in refobj.references:
        benchmark_stats["numnonexcludedbases"][ref] = refobj.get_reference_length(ref)
    if bedobjects["allexcludedregions"] is not None:
        for interval in bedobjects["allexcludedregions"]:
            ref = interval.chrom
            if ref in benchmark_stats["numnonexcludedbases"].keys():
                benchmark_stats["numnonexcludedbases"][ref] = benchmark_stats["numnonexcludedbases"][ref] - len(interval)
    
    aligndict = {}
    for align in aligndata:
        refentry = align["target"]
        if refentry not in aligndict:
            aligndict[refentry] = [align]
        else:
            aligndict[refentry].append(align)

    # assess each benchmark entry, one by one
    for refentry in sorted(aligndict.keys()):
        logger.debug("Entry " + refentry + " start")
        refnelength = benchmark_stats["numnonexcludedbases"][refentry]
        refalignclusters = []
        # sort alignments from longest (along the benchmark) to shortest:
        numaligns = len(aligndict[refentry])
        logger.debug("Sorting " + str(numaligns) + " " + refentry + " aligns")
        aligndict[refentry].sort(reverse=True, key=lambda align: align["targetalignlength"])
        # calculate slope as query diff over ref diff, and cluster alignments within the same query/target band:
        logger.debug("Adding " + str(numaligns) + " " + refentry + " aligns to clusters")
        for refalign in aligndict[refentry]:
            alignrefstart = refalign['targetstart']
            alignrefend = refalign['targetend']
            alignquery = refalign['query']
            alignquerystart = refalign['querystart']
            alignqueryend = refalign['queryend']
            add_align_to_clusters(refalign, refalignclusters, maxdistance)

        # split clusters that are separated along the target by more than maxdistance:
        logger.debug("Splitting " + refentry + " clusters")
        refalignclusters = split_disjoint_clusters(refalignclusters, maxdistance)

        logger.debug("Forming " + refentry + " cluster bed string")
        for cluster in sorted(refalignclusters, key=lambda c: c["aligns"][0]["targetstart"]):
            clusterquery = cluster["query"] 
            clusterslope = cluster["slope"]
            clusterintercept = cluster["intercept"]
            clusterbedstring = ""
            for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
                clusterbedstring = clusterbedstring + refentry + "\t" + str(align["targetstart"]) + "\t" + str(align["targetend"]) + "\t" + clusterquery + "_" + str(align["querystart"]) + "_" + str(align["queryend"]) + "\n"

            bedtool = pybedtools.BedTool(clusterbedstring, from_string = True)
            mergedbedtool = bedtoolslib.mergeintervals(bedtool)
            clusterbases = bedtoolslib.bedsum(mergedbedtool)
            if bedobjects["allexcludedregions"] is not None:
                nonexcludedbedtool = bedtoolslib.subtractintervals(mergedbedtool, bedobjects["allexcludedregions"])
            else:
                nonexcludedbedtool = mergedbedtool
            cluster["nonexcludedcoveredbases"] = bedtoolslib.bedsum(nonexcludedbedtool)
            logger.debug("Cluster on " + clusterquery + " has " + str(clusterbases) + " non-redundant bases, " + str(cluster["nonexcludedcoveredbases"]) + " of which are not excluded")

        logger.debug("See how many " + refentry + " clusters are needed to cover 95% of ref")
        # calculate how many clusters needed to cover 95% of ref:
        totalnonexcludedcovered = 0
        clustercount = 0
        lca95 = None
        nca95 = None
        clusterno = 1
        refentrybedstring = ""
        for cluster in sorted(refalignclusters, key=lambda c:c["nonexcludedcoveredbases"], reverse = True):
            clusterquery = cluster["query"] 
            for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
                if lca95 is None:
                    clustername = "Cluster" + str(clusterno)
                else:
                    clustername = "SmallCluster" + str(clusterno)
                refentrybedstring = refentrybedstring + refentry + "\t" + str(align["targetstart"]) + "\t" + str(align["targetend"]) + "\t" + clusterquery + "_" + str(align["querystart"]) + "_" + str(align["queryend"]) + "_" + clustername + "\n"
            clusterno = clusterno + 1

            nebases = cluster["nonexcludedcoveredbases"]
            logger.debug("Cluster " + clustername + " has " + str(nebases) + " nonexcluded covered bases")

            if lca95 is None:
                totalnonexcludedcovered = totalnonexcludedcovered + cluster["nonexcludedcoveredbases"]
            clustercount = clustercount + 1
            if lca95 is None and totalnonexcludedcovered > 0.95*refnelength:
                lca95 = clustercount
                nca95 = cluster["nonexcludedcoveredbases"]

        logger.debug("Saving BED file for " + refentry + " clusters")
        refentrybedtool = pybedtools.BedTool(refentrybedstring, from_string = True)
        sortedrefentrybedtool = refentrybedtool.sort()
        sortedrefentrybedtool.saveas(alignplotprefix + "." + refentry + ".clusters.bed")

        # Note: lca95 will be "None" for reference entries not able to be covered
        benchmark_stats["clustercoverage"][refentry] = {"lca95":lca95, "nonexcludedcovered":totalnonexcludedcovered}
        refalignclusters.sort(key=lambda c: c["aligns"][0]["targetstart"])
        benchmark_stats["alignclusters"][refentry] = refalignclusters
    
    return benchmark_stats

def left_right_soft_clip(align)->list:
    leftsoftclip = 0
    rightsoftclip = 0
    
    if align.cigartuples[0][0]==4:
        leftsoftclip = align.cigartuples[0][1]
    if align.cigartuples[-1][0]==4:
        rightsoftclip = align.cigartuples[-1][1]

    return [leftsoftclip, rightsoftclip]

def left_right_hard_clip(align)->list:
    lefthardclip = 0
    righthardclip = 0
    
    if align.cigartuples[0][0]==5:
        lefthardclip = align.cigartuples[0][1]
    if align.cigartuples[-1][0]==5:
        righthardclip = align.cigartuples[-1][1]

    return [lefthardclip, righthardclip]

# Routine to retrieve alignment ref coords and cigar ops for subalignment limits (specified by offsets from query_alignment_start)
def retrieve_refcoords_and_cigars_from_querycoords(align, desiredquerystart:int, desiredqueryend:int):

    # general stats about this alignment:
    query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
    #print("Align " + query + ":" + str(querystart) + "-" + str(queryend) + " looking for positions " + str(desiredquerystart) + " and " + str(desiredqueryend))
    [left_softclip, right_softclip] = left_right_soft_clip(align)
    [left_hardclip, right_hardclip] = left_right_hard_clip(align)
    alignops = align.cigartuples

    # first position will begin at left-most ref/query base of the alignment (regardless of strand)
    refcurrentoffset = 0
    querycurrentoffset = 0
    currentalignopindex = 0

    desiredrefstart = None
    desiredrefend = None
    desiredcigarops = []
    desiredfirstopindex = None
    desiredlastopindex = None
    while refcurrentoffset <= refend-refstart and currentalignopindex < len(alignops): # traverse the alignment operator by operator
        alignop = alignops[currentalignopindex]
        op = alignop[0]
        oplength = alignop[1]

        currentalignopindex = currentalignopindex + 1

        if op in [0, 7, 8]: # substitution
            matches_to_add = oplength
            if desiredquerystart >= querycurrentoffset and desiredquerystart < querycurrentoffset + oplength:
                offsetwithinblock = desiredquerystart - querycurrentoffset
                desiredrefstart = refcurrentoffset + refstart + offsetwithinblock
                matches_to_add = matches_to_add - offsetwithinblock 
            if desiredqueryend >= querycurrentoffset and desiredqueryend < querycurrentoffset + oplength:
                offsetwithinblock = desiredqueryend - querycurrentoffset
                desiredrefend = refcurrentoffset + refstart + offsetwithinblock
                matches_to_add = matches_to_add - oplength + offsetwithinblock 
            if desiredrefstart is not None:
                desiredcigarops.append([op, matches_to_add])
            if desiredrefend is not None:
                return desiredquerystart, desiredqueryend, ref, desiredrefstart, desiredrefend, desiredcigarops
            refcurrentoffset = refcurrentoffset + oplength
            querycurrentoffset = querycurrentoffset + oplength
        # advance current positions: cases where reference coord advances (DN):
        elif op in [2, 3]:
            if desiredrefstart is not None:
                desiredcigarops.append([op, oplength])
            refcurrentoffset = refcurrentoffset + oplength
        # cases where query coord advances (I)
        # BUG: need to account here for cases where desired query endpoints are reached within an insertion
        elif op == 1:
            insertions_to_add = oplength
            if desiredquerystart >= querycurrentoffset and desiredquerystart < querycurrentoffset + oplength:
                desiredrefstart = refcurrentoffset + refstart
                insertions_to_add = insertions_to_add + querycurrentoffset - desiredquerystart
            if desiredqueryend >= querycurrentoffset and desiredqueryend < querycurrentoffset + oplength:
                desiredrefend = refcurrentoffset + refstart
                insertions_to_add = desiredqueryend - querycurrentoffset

                desiredcigarops.append([op, insertions_to_add])

                return desiredquerystart, desiredqueryend, ref, desiredrefstart, desiredrefend, desiredcigarops

            if desiredrefstart is not None:
                desiredcigarops.append([op, insertions_to_add])
            querycurrentoffset = querycurrentoffset + oplength

        if desiredrefstart is not None and desiredrefend is not None:
            return desiredquerystart, desiredqueryend, ref, desiredrefstart, desiredrefend, desiredcigarops

    return desiredquerystart, desiredqueryend, ref, desiredrefstart, desiredrefend, desiredcigarops

# do a quick comparison of two alignments to see if one is contained within the other, or if they are dovetail to the left, right, or both
# function returns either: ref start/end within the first alignment where the second alignment has agreeing endpoints, but -1 for 
# start or end if there is no agreement at that end
def compare_alignments(align1, align2, switched=False):

    # coordinates are all one-based, with start at beginning of *original* sequence (not left end of the alignment)
    query1, query1start, query1end, ref1, ref1start, ref1end, strand1 = retrieve_align_data(align1)
    query2, query2start, query2end, ref2, ref2start, ref2end, strand2 = retrieve_align_data(align2)

    logger.debug("Comparing " + ref1 + ":" + str(ref1start) + "-" + str(ref1end) + " to " + query2 + ":" + str(query2start) + "-" + str(query2end) + "/" + strand1)

    if strand1 != strand2:
        return [-1, -1]

    # for soft clipped alignments, query positions include soft clipped bases (so are larger by the amount of the soft clipping)
    align1pairs = align1.get_aligned_pairs(matches_only=True)
    align2pairs = align2.get_aligned_pairs(matches_only=True)

    [leftsoftclip1, rightsoftclip1] = left_right_soft_clip(align1)
    [leftsoftclip2, rightsoftclip2] = left_right_soft_clip(align2)

    # for align1, ref positions in array are the second tuple element + ref1start, query pos are first tuple element + query1start
    align1index = 0
    align2index = 0
    alignoverlaprefstart = None
    alignoverlapquerystart = None
    lastoverlaprefstart = None
    lastoverlapquerystart = None
    while align1index < len(align1pairs) and align2index < len(align2pairs):
        ref1pos = align1pairs[align1index][1] + 1
        if strand1 == 'F':
            query1pos = align1pairs[align1index][0] + query1start
            if leftsoftclip1 > 0:
                query1pos = query1pos - leftsoftclip1
        else:
            query1pos = align1pairs[align1index][0] + query1start
            if rightsoftclip1 > 0:
                query1pos = query1pos - rightsoftclip1
    
        if not switched:
            ref2pos = align2pairs[align2index][1] + 1
            if strand2 == 'F':
                query2pos = align2pairs[align2index][0] + query2start
                if leftsoftclip2 > 0:
                    query2pos = query2pos - leftsoftclip2
            else:
                query2pos = query2start + align2pairs[align2index][0]
                if rightsoftclip2 > 0:
                    query2pos = query2pos - rightsoftclip2

        else: # need to switch ref/query
            if strand2 == 'F':
                ref2pos = align2pairs[align2index][0] + query2start
                query2pos = align2pairs[align2index][1] + 1
                if leftsoftclip2 > 0:
                    ref2pos = ref2pos - leftsoftclip2
            else:
                ref2pos = align2pairs[align2index][0] + query2start
                #query2pos = ref2end - align2pairs[align2index][1]
                query2pos = align2pairs[align2index][1] + 1
                if rightsoftclip2 > 0:
                    ref2pos = ref2pos - rightsoftclip2

        # based on ref positions, should we advance one index, the other, or both?
        #if ref1start==21156327 and ref1end==34675573:
            #print(str(ref1pos) + "\t" + str(ref2pos) + "\t" + str(query1pos) + "\t" + str(query2pos) + "\t" + str(align1index) + "\t" + str(align2index))
        if ref1pos < ref2pos:
            align1index = align1index + 1
        elif ref1pos > ref2pos:
            align2index = align2index + 1
        elif ref1pos == ref2pos:
            if query1pos == query2pos:
                if alignoverlaprefstart is None:
                    alignoverlaprefstart = ref1pos
                    alignoverlapquerystart = query1pos
                lastoverlaprefstart = ref1pos
                lastoverlapquerystart = query1pos
            align1index = align1index + 1
            align2index = align2index + 1

    return [alignoverlaprefstart, lastoverlaprefstart]

def add_align_to_clusters(align:dict, alignclusters:list, maxdistance:int):

    alignstart = align['targetstart']
    alignend = align['targetend']
    alignquery = align['query']
    alignquerystart = align['querystart']
    alignqueryend = align['queryend']
    alignslope = (alignend - alignstart)/(alignqueryend-alignquerystart)
    alignintercept = alignstart - int(alignslope * alignquerystart)

    # try to assign this align to a pre-existing cluster of aligns:
    assigned = False
    for cluster in alignclusters:
        if cluster["query"] != alignquery:
            continue
        clusterquery = cluster["query"]
        clusterslope = cluster["slope"]
        clusterintercept = cluster["intercept"]
        predstart = clusterintercept + clusterslope * alignquerystart
        if abs(predstart - alignstart) <= maxdistance:
            cluster["aligns"].append(align)
            assigned = True
            break

    # create a new cluster if none were appropriate
    if not assigned:
        alignclusters.append({'query':alignquery, 'slope':alignslope, 'intercept':alignintercept, 'aligns':[align]})

    return 0

def split_disjoint_clusters(refalignclusters:list, maxdistance:int)->list:

    spinoffclusters = []
    for cluster in refalignclusters:
        if len(cluster["aligns"])==1:
            continue

        runningalignlist = []
        maxposition = 0
        for align in sorted(cluster["aligns"], key=lambda a: (a["targetstart"], a["targetend"])):
            # do we have too big a gap?
            if len(runningalignlist) > 0 and align["targetstart"] - maxposition > maxdistance:
                spinoffclusters.append({'query':cluster["query"], 'slope':cluster["slope"], 'intercept':cluster["intercept"], 'aligns':runningalignlist})
                runningalignlist = []
                maxposition = 0
                logger.debug("Splitting align cluster that ends at " + cluster["query"] + ":" + str(maxposition) + " for new align starting at " + str(align["targetstart"]))
            runningalignlist.append(align)
            maxposition = max(maxposition, align["targetend"])
        cluster["aligns"] = runningalignlist

    return refalignclusters + spinoffclusters

def trim_bamfile_to_intervals(bamfile, intervals, outputbam, headerbam, args, sort=True, index=True):

    [aligns, alignedintervals] = index_aligns_by_boundaries(bamfile, args)
    subaligns = find_phaseblock_subaligns(intervals, alignedintervals, aligns)
    write_aligns_to_bamfile(outputbam, subaligns, headerbam=headerbam, sort=sort, index=index)

def index_aligns_by_boundaries(bamfile, args):

    logger.info("Indexing aligns in " + bamfile)
    alignobj = pysam.AlignmentFile(bamfile, "rb")
    aligndict = {}
    alignbedstring = ""
    for align in alignobj.fetch():
        if align.is_secondary:
            continue
        if align.reference_length >= args.minalignlength:
            query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(align)
            alignname = query + "_" + str(querystart) + "_" + str(queryend)
            alignbedstring = alignbedstring + query + "\t" + str(querystart) + "\t" + str(queryend) + "\t" + str(alignname) + "\n"
            aligndict[alignname] = align

    alignbedtool = pybedtools.BedTool(alignbedstring, from_string = True)

    return [aligndict, alignbedtool]

# Trim alignments to the boundaries of phase blocks:
def find_phaseblock_subaligns(phaseblockints, alignedintervals, aligndict):

    # For each alignment, create a list of intersected intervals representing phaseblock regions within the alignment
    # Intervals are in the direction of the query, even if the alignment is on the reverse strand
    phasedaligndict = {}
    phaseblockalignints = bedtoolslib.intersectintervals(phaseblockints, alignedintervals, wo=True)
    numphaseblockaligns = len(phaseblockalignints)
    logger.debug("There are " + str(numphaseblockaligns) + " aligned intervals within phase block intervals")
    for blockalign in phaseblockalignints:
        alignname = blockalign[12]
        if alignname not in phasedaligndict.keys():
            phasedaligndict[alignname] = []
        phasedaligndict[alignname].append(blockalign)

    # Now for each aligned region, sort the intersecting intervals from low to high coordinate:
    sortedaligns = sorted(phasedaligndict.keys())
    subalignlist = []
    for alignname in sortedaligns:
        if alignname not in aligndict.keys():
            logger.debug("Skipping alignment " + alignname + "--not in align dictionary!")
            print("Skipping alignment " + alignname + "--not in align dictionary!")
            continue
        alignobj = aligndict[alignname]
        query, querystart, queryend, ref, refstart, refend, strand = retrieve_align_data(alignobj)
        phaseblockintersects = phasedaligndict[alignname]
        subaligninfo = []
        segnumber = 1
        logger.debug("Processing align " + alignname)
        for phaseblock in sorted(phaseblockintersects, key=lambda x:(x.start, x.stop)):
            pbstart = phaseblock.start
            pbend = phaseblock.end
            alignstart = int(phaseblock[10])
            alignend = int(phaseblock[11])
            intersectlength = int(phaseblock[13])
            intersectstart = max(pbstart, alignstart)
            intersectend = min(pbend, alignend)
            logger.debug("Align " + alignname + " intersects phase block " + str(pbstart) + "-" + str(pbend) + " from " + str(intersectstart) + " to " + str(intersectend))
            
            if (intersectend - intersectstart) != intersectlength:
                print(str(intersectlength) + " is not equal to " + str(intersectend) + "-" + str(intersectstart))
            # Need to calculate the correct aligned query start/end for the routine "create_subalignobjects":
            if strand == 'F':
                querystartoffset = intersectstart - querystart
                queryendoffset = intersectend - querystart
            else: # this should be corrected:
                querystartoffset = queryend - intersectend
                queryendoffset = queryend - intersectstart

            desiredquerystart, desiredqueryend, desiredref, desiredrefstart, desiredrefend, desiredcigarops = retrieve_refcoords_and_cigars_from_querycoords(alignobj, querystartoffset, queryendoffset)

            desiredcigarqueryconsumed = count_consumed_query(desiredcigarops)
            if desiredcigarqueryconsumed != intersectlength:
               print("Align " + alignname + " intersects phase block " + str(pbstart) + "-" + str(pbend) + " from " + str(intersectstart) + " to " + str(intersectend))
               print("Found refcoords " + desiredref + ":" + str(desiredrefstart) + "-" + str(desiredrefend))
               print("Cigar ops consume " + str(desiredcigarqueryconsumed) + " bases and intersect is " + str(intersectlength))
               continue

            if desiredrefstart is None or desiredrefend is None:
                print("Found refcoords " + desiredref + ":" + str(desiredrefstart) + "-" + str(desiredrefend))
                print("Unable to find reference start/end for query " + str(querystartoffset) + "-" + str(queryendoffset) + " in align " + alignname)
                continue
            subaligninfo.append({'alignedquerystart':querystartoffset, 'alignedqueryend':queryendoffset, 'alignedrefstart':desiredrefstart, 'alignedrefend':desiredrefend, 'cigarops':desiredcigarops, 'subalignlength':intersectlength, 'segnum':segnumber})
            segnumber = segnumber + 1

        [left_hardclip, right_hardclip] = left_right_hard_clip(alignobj)
        hardcliplongest = False
        if alignobj.is_supplementary or (left_hardclip > 0) or (right_hardclip > 0):
            hardcliplongest = True
        subalignobjs = create_subalignobjects(alignobj, subaligninfo, hardcliplongest)
        subalignlist = subalignlist + subalignobjs

    return subalignlist

def write_aligns_to_bamfile(bamfilename:str, aligns:list, headerbam:str="test_bam_for_header.bam", sort=True, index=True):

    bamobj = pysam.AlignmentFile(headerbam, "r")

    with pysam.AlignmentFile(bamfilename, "wb", header=bamobj.header) as sbfh:
        for align in aligns:
            sbfh.write(align)

    if sort:
        sortedbamfilename = bamfilename.replace(".bam", ".sort.bam")
        pysam.sort("-o", sortedbamfilename, bamfilename)
        if index:
            pysam.index(sortedbamfilename)
        return sortedbamfilename
    else:
        return bamfilename

def write_aligns_to_samfile(samfilename:str, aligns:list, headerbam:str="test_bam_for_header.bam"):

    bamobj = pysam.AlignmentFile(headerbam, "r")

    print("Writing alignments to " + samfilename)
    with pysam.AlignmentFile(samfilename, "w", header=bamobj.header) as sbfh:
        for align in aligns:
            sbfh.write(align)


def merge_trimmed_bamfiles(mattrimmedbamfile:str, pattrimmedbamfile:str, benchdiploidheaderfile:str, outputfiles:dict, sort=True):

    mergedtrimmedsamfile = outputfiles["trimmedphasedalignprefix"] + ".merge.dipheader.sam"
    headerreturnval = os.system("cat " + benchdiploidheaderfile + " > " + mergedtrimmedsamfile)
    logger.debug("Copying " + benchdiploidheaderfile + " to " + mergedtrimmedsamfile + "had return value " + str(headerreturnval))

    mattrimmedsamreturnval = os.system("samtools view " + mattrimmedbamfile + " >> " + mergedtrimmedsamfile)
    logger.debug("Writing maternal trimmed sam alignments to " + mergedtrimmedsamfile + " had return value " + str(mattrimmedsamreturnval))
    pattrimmedsamreturnval = os.system("samtools view " + pattrimmedbamfile + " >> " + mergedtrimmedsamfile)
    logger.debug("Writing paternal trimmed sam alignments to " + mergedtrimmedsamfile + " had return value " + str(pattrimmedsamreturnval))
    mergedtrimmedbamfile = outputfiles["trimmedphasedalignprefix"] + ".merge.bam"
    samtobamreturnval = os.system("samtools view -bS " + mergedtrimmedsamfile + " -o " + mergedtrimmedbamfile)
    logger.debug("Converting merged trimmed sam alignments to bam had return value " + str(samtobamreturnval))

    if sort:
        trimmedphasedsortbam = outputfiles["trimmedphasedalignprefix"] + ".merge.sort.bam"
        trimmedbamsortreturnval = os.system("samtools sort " + mergedtrimmedbamfile + " > " + trimmedphasedsortbam)
        logger.debug("Sorting merged trimmed bam alignments had return value " + str(trimmedbamsortreturnval))
        trimmedbamindexval = os.system("samtools index " + trimmedphasedsortbam)
        logger.debug("Indexing trimmed bam file had return value " + str(trimmedbamindexval))
        if trimmedbamsortreturnval==0 and trimmedbamindexval==0:
            os.remove(mergedtrimmedsamfile)
            return ""
        return trimmedphasedsortbam

    return mergedtrimmedbamfile
