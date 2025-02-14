import sys
import logging

logger = logging.getLogger(__name__)

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
            if current_align is not None and refentry == current_align["target"]:
                if refstart < current_align["targetend"]:
                    if query == current_align["query"]:
                        sfh.write(refentry + "\t" + str(refstart) + "\t" + str(current_align["targetend"]) + "\tSameContigInsertion\n")
                    else:
                        sfh.write(refentry + "\t" + str(refstart) + "\t" + str(current_align["targetend"]) + "\tBetweenContigInsertion\n")
                elif refstart >= current_align["targetend"]:
                    if query == current_align["query"]:
                        sfh.write(refentry + "\t" + str(align["targetend"]) + "\t" + str(refstart) + "\tSameContigDeletion\n")
                    else:
                        sfh.write(refentry + "\t" + str(align["targetend"]) + "\t" + str(refstart) + "\tBetweenContigDeletion\n")
            current_align = align

    return 0

