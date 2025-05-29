import sys
import logging

logger = logging.getLogger(__name__)

# reminders: in aligndata, (1) all coordinates are 1-based, (2) strand is "+" or "-", (3) querystart is the query's lower coordinate,
# so doesn't correspond to targetstart if alignment is on the reverse strand

def write_structural_errors(aligndata:list, refobj, queryobj, outputdict, bmstats, args)->str:

    aligndict = {}
    current_align = None
    with open(outputdict["structvariantbed"], "w") as sfh:
        for align in sorted(aligndata, key=lambda a: (a["target"], a["targetstart"], a["targetend"])):
            refentry = align["target"]
            if refentry not in aligndict:
                aligndict[refentry] = [align]
            else:
                aligndict[refentry].append(align)
            query = align["query"]
            refstart = align["targetstart"]
            refend = align["targetend"]
            querystart = align["querystart"]
            queryend = align["queryend"]
            strand = align["strand"]
            if current_align is not None:
                refdiff = refstart - current_align["targetend"]
                if refentry == current_align["target"] and query == current_align["query"] and strand == current_align["strand"]:
                    if strand == "+":
                        querydiff = querystart - current_align["queryend"]
                        query1 = current_align["queryend"]
                        query2 = querystart
                    else:
                        querydiff = queryend - current_align["querystart"]
                        query1 = querystart
                        query2 = current_align["queryend"]
   
                    netdiff = querydiff - refdiff
                    if refdiff < querydiff: # refdiff less than querydiff (insertion), netshift positive
                        if refdiff > 0:
                            sfh.write(refentry + "\t" + str(current_align["targetend"] - 1) + "\t" + str(refstart) + "\tSameContigInsertion\t" + query + "\t" + str(query1) + "\t" + str(query2) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n")
                        else:
                            sfh.write(refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tSameContigInsertion\t" + query + "\t" + str(query1) + "\t" + str(query2) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n")
                    else: # refdiff greater than than querydiff (deletion), netshift negative
                        if refdiff > 0:
                            sfh.write(refentry + "\t" + str(current_align["targetend"] - 1) + "\t" + str(refstart) + "\tSameContigDeletion\t" + query + "\t" + str(query1) + "\t" + str(query2) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n")
                        else:
                            sfh.write(refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tSameContigDeletion\t" + query + "\t" + str(query1) + "\t" + str(query2) + "\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\t" + str(netdiff) + "\t" + strand + "\n")
    
                elif refentry == current_align["target"]: # strand switch or new contig:
                    queryentries = query + "/" + current_align["query"]
                    strands = strand + "/" + current_align["strand"]
                    if refdiff > 0:
                        sfh.write(refentry + "\t" + str(current_align["targetend"] - 1) + "\t" + str(refstart) + "\tBetweenContigDeletion\t" + queryentries + "\t.\t.\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\tNA\t" + strands + "\n")
                    else:
                        sfh.write(refentry + "\t" + str(refstart - 1) + "\t" + str(current_align["targetend"]) + "\tBetweenContigInsertion\t" + queryentries + "\t.\t.\t" + str(current_align["targetend"]) + "\t" + str(refstart) + "\tNA\t" +strands + "\n")
            current_align = align

    return 0

