import re
import pysam
import pybedtools
import logging
from collections import namedtuple
from GQC import seqparse
from GQC import phasing

logger = logging.getLogger(__name__)

# This routine, ported from mummer's delta-filter C++ routines "flagQLIS" and "flagRLIS", finds the longest
# increasing sequence of alignments along a target- or query-sorted list of alignments, using the
# method's "sorttype" argument to determine whether to do QLIS (sorttype="query") or RLIS (sorttype=
# "target").
#
# Inputs: The list of alignments passed to the routine should contain a dictionary for each align
# with fields: target, targetstart, targetend, targetlength, query, querystart, queryend, querylength
# For reverse strand alignments, querystart can be larger than queryend, but targetstart should not 
# be larger than targetend. All coordinates are assumed to be 1-based.
#
# Return value: The function returns a list with the same alignments passed as input, except the list 
# doesn't contain any of the aligns that are part of the QLIS or RLIS, and it is sorted in order of
# increasing query or reference lower-coordinate.
#

def filter_aligns(alignlist:list, sorttype="target", maxoverlap=0.95, mummeralgorithm=True):

    # define the fields that will be used in sorting and determination of length, etc.
    if sorttype == "target":
        sortfield = "targetstart"
        lowfield = "targetstart"
        highfield = "targetend"
    else:
        sortfield = "querylow"
        lowfield = "querylow"
        highfield = "queryhigh"
   
    for align in alignlist:
        # set query high/low
        if align["querystart"] < align["queryend"]:
            align["querylow"] = align["querystart"]
            align["queryhigh"] = align["queryend"]
        else:
            align["querylow"] = align["queryend"]
            align["queryhigh"] = align["querystart"]
    
    # make a hash grouping the alignments by target or query entry
    alignbyentrydict = create_align_by_entry_dict(alignlist, entrytype=sorttype)

    # create a hash with keys for each alignment that passes the LIS filter
    filteredaligns = []
    for entry in alignbyentrydict:
        alignlist = sorted(alignbyentrydict[entry], key=lambda align: align[sortfield])
        numaligns = len(alignlist)
        lis = []

        if mummeralgorithm:
            for i in range(numaligns):
                lis.append({"used":False})
            allbest = []
    
            while True:
                for i in range(numaligns):
                    if lis[i]["used"]:
                        continue
                    if "identity" not in alignlist[i]:
                        alignlist[i]["identity"] = 1.0
                    lis[i]["align"] = alignlist[i]
                    lis[i]["score"] = (alignlist[i][highfield] - alignlist[i][lowfield] + 1)*alignlist[i]["identity"]*alignlist[i]["identity"]
                    lis[i]["from"] = -1
                    lis[i]["diff"] = 0
    
                    for j in range(i):
                        if lis[j]["used"] or (lis[j]["from"] >=0 and lis[lis[j]["from"]]["align"][highfield] >= lis[i]["align"][lowfield]):
                            continue
                        leni = lis[i]["align"][highfield] - lis[i]["align"][lowfield] + 1
                        lenj = lis[j]["align"][highfield] - lis[j]["align"][lowfield] + 1
                        olap = lis[j]["align"][highfield] - lis[i]["align"][lowfield] + 1
                        if olap < 0:
                            olap = 0
                        diff = lis[j]["diff"]
                        if lis[j]["align"]["targetstart"] < lis[i]["align"]["targetstart"]:
                            diff = diff + lis[j]["align"]["targetend"] - lis[i]["align"]["targetstart"]
                        else:
                            diff = diff + lis[i]["align"]["targetend"] - lis[j]["align"]["targetstart"]
                        if lis[j]["align"]["querylow"] < lis[i]["align"]["querylow"]:
                            diff = diff + lis[j]["align"]["queryhigh"] - lis[i]["align"]["querylow"]
                        else:
                            diff = diff + lis[i]["align"]["queryhigh"] - lis[j]["align"]["querylow"]
    
                        olapfraction = max(olap/leni, olap/lenj)
                        if olapfraction <= maxoverlap:
                            score = lis[j]["score"] + (leni - olap)*alignlist[i]["identity"]*alignlist[i]["identity"]
                        else:
                            score = -1.0
    
                        if score > lis[i]["score"] or (score == lis[i]["score"] and diff < lis[i]["diff"]):
                            lis[i]["from"] = j
                            lis[i]["score"] = score
                            lis[i]["diff"] = diff
                if updatebest(lis, numaligns, allbest):
                    break
    
            numbests = len(allbest)
            eqc = 0
            while eqc < numbests:
                if lis[allbest[eqc]]["diff"] != lis[allbest[0]]["diff"]:
                    break
                eqc = eqc + 1
            # need to have this pick a random number between 0 and eqc instead of using 0
            bestpick = 0
    
            i = allbest[bestpick]
            while i != -1:
                filteredaligns.insert(0, lis[i]["align"])
                i = lis[i]["from"]
        else:
            return -1
        
    return filteredaligns
        #foreach my $rh_entry_pair (@{$ra_entry_pairs}) {
            #my $refentry = $rh_entry_pair->{ref_entry};
            #my $queryentry = $rh_entry_pair->{query_entry};
            #print STDERR "Flagging pair $refentry $queryentry\n";
            #my $ra_aligns = $rh_entry_pair->{aligns};
            #foreach my $rh_align (@{$ra_aligns}) {
                #$rh_align->{query_redundant} = 1;
                #$rh_align->{ref_redundant} = 1;
            #}
    #
            ## find and mark query_redundant:
            #my @sorted_aligns = sort {$a->{query_low} <=> $b->{query_low}} @{$ra_aligns};
    #
            ## mark alignments that are non-redundant by query coord:
            #my @query_starts = map {$_->{query_start}} @{$ra_aligns};
            #my @query_ends = map {$_->{query_end}} @{$ra_aligns};
    #
            ## sorted list of all endpoints:
            #my $no_query_ends = @query_ends;
            #my $no_query_starts = @query_starts;
            #my $expected_total = $no_query_ends + $no_query_starts;
            #my @all_query_ends = @query_starts;
            #push @all_query_ends, @query_ends;
    #
            #@all_query_ends = sort {$a <=> $b} @all_query_ends;
    #
            #my $actual_total = @all_query_ends;
            #print STDERR "Expected total $expected_total actual total $actual_total\n";
    #
            #my @sweep = ();
            #my $last_end = 0;
            #for (my $i=0; $i<=$#all_query_ends; $i++) {
                #my $no_sweep = @sweep;
                #print STDERR "At endpoint $i with $no_sweep in sweep\n";
                #my $this_end = $all_query_ends[$i];
                #next if ($this_end == $last_end); # don't repeat end values
                #$last_end = $this_end; # new end to process
                #while( (@sorted_aligns) && ($sorted_aligns[0]->{query_low}==$this_end )) {
                    #push @sweep, shift @sorted_aligns;
                #}
    #
                #my @temp_sweep = ();
                #my $max_score = 0;
                #foreach my $rh_align (@sweep) {
                    #if ($rh_align->{query_high} != $this_end) {
                        #push @temp_sweep, $rh_align;
                        #if (($rh_align->{query_high} - $rh_align->{query_low}) > $max_score) {
                            #$max_score = $rh_align->{query_high} - $rh_align->{query_low};
                        #}
                    #}
                #}
                #@sweep = @temp_sweep; # keeps only segs that don't end at this end
                #foreach my $rh_align (@sweep) { # mark max score aligns as "good"--scores should have been stored in hash
                    #if (($rh_align->{query_high} - $rh_align->{query_low}) == $max_score) {
                        #$rh_align->{query_redundant} = 0;
                    #}
                #}
            #}
    #
            ## find and mark ref_redundant:
            #@sorted_aligns = sort {$a->{ref_start} <=> $b->{ref_start}} @{$ra_aligns};
    #
            ## mark alignments that are non-redundant by ref coord:
            #my @ref_starts = map {$_->{ref_start}} @{$ra_aligns};
            #my @ref_ends = map {$_->{ref_end}} @{$ra_aligns};
    #
            ## sorted list of all endpoints:
            #my @all_ref_ends = @ref_starts;
            #push @all_ref_ends, @ref_ends;
            #@all_ref_ends = sort {$a <=> $b} @all_ref_ends;
    #
            #@sweep = ();
            #$last_end = 0;
            #for (my $i=0; $i<=$#all_ref_ends; $i++) {
                #my $no_sweep = @sweep;
                #print STDERR "At endpoint $i with $no_sweep in sweep\n";
                #my $this_end = $all_ref_ends[$i];
                #next if ($this_end == $last_end); # don't repeat end values
                #$last_end = $this_end;
                #while( (@sorted_aligns) && ($sorted_aligns[0]->{ref_start}==$this_end )) {
                    #push @sweep, shift @sorted_aligns;
                #}
    #
                #my @temp_sweep = ();
                #my $max_score = 0;
                #foreach my $rh_align (@sweep) {
                    #if ($rh_align->{ref_end} != $this_end) {
                        #push @temp_sweep, $rh_align;
                        #if (($rh_align->{ref_end} - $rh_align->{ref_start}) > $max_score) {
                            #$max_score = $rh_align->{ref_end} - $rh_align->{ref_start};
                        #}
                    #}
                #}
                #@sweep = @temp_sweep; # keeps only segs that don't end at this end
                #foreach my $rh_align (@sweep) { # mark max score aligns as "good"--scores should have been stored in hash
                    #if (($rh_align->{ref_end} - $rh_align->{ref_start}) == $max_score) {
                        #$rh_align->{ref_redundant} = 0;
                    #}
                #}
            #}
        #}
    #}


def create_align_by_entry_dict(alignlist:list, entrytype:str):

    alignbyentrydict = {}
    for align in alignlist:
        entry = align[entrytype]
        if entry not in alignbyentrydict:
            alignbyentrydict[entry] = [align]
        else:
            alignbyentrydict[entry].append(align)

    return alignbyentrydict

def updatebest(lis:list, n:int, allbest:list):
    if n==0:
        return False
    for best in range(n):
        if not lis[best]["used"]:
            break
    for i in range(best+1, n, 1):
        if not lis[i]["used"] and (lis[i]["score"] > lis[best]["score"] or (lis[i]["score"]==lis[best]["score"] and lis[i]["diff"] < lis[best]["diff"]) ):
            best = i
    if len(allbest) > 0 and lis[allbest[0]]["score"] > lis[best]["score"]:
        return False

    allbest.append(best)
    i = best
    while i != -1:
        lis[i]["used"] = True
        i = lis[i]["from"]

    return True

# This routine, ported from mummer's show-diff C++ routine "PrintDiff", steps through reference
# ordered alignments (which have been previously filtered to be RLIS), and analyses the 
# alignment breaks to determine what type they are.
#
# Inputs: The list of alignments passed to the routine should contain a dictionary for each align
# with fields: target, targetstart, targetend, targetlength, query, querystart, queryend, querylength,
# and strand, which should be "F" or "R". All coordinates are assumed to be 1-based.
#
# Return value: The function returns a list with the same alignments passed as input, except the list 
# doesn't contain any of the aligns that are part of the QLIS or RLIS, and it is sorted in order of
# increasing query or reference lower-coordinate.
#

def find_diffs(alignlist:list):
    
    alignbyentrydict = create_align_by_entry_dict(alignlist, entrytype="target")

    inversionends = []
    lisjumps = []
    for target in alignbyentrydict.keys():
        entryaligns = alignbyentrydict[target]
        prevalign = None
        for align in entryaligns:
            # no breaks possible in first align, so just record as prevalign and go to second:
            if not prevalign:
                prevalign = align
                continue
            # break analysis
            gapR = align["targetstart"] - prevalign["targetend"] - 1
            prevquery = prevalign["query"]
            alignquery = align["query"]

            # breaks between aligns to the same scaffold:
            if prevquery == alignquery:
                if prevalign["strand"] != align["strand"]:
                    inversionends.append({"target":target, "query":query, "prevstrand":prevalign["strand"], "alignstrand":align["strand"], "prevtargetend":prevalign["targetend"], "aligntargetstart":align["targetstart"]})
                else:
                    lisjumps.append({"target":target, "query":query, "prevstrand":prevalign["strand"], "alignstrand":align["strand"], "prevtargetend":prevalign["targetend"], "aligntargetstart":align["targetstart"], "prevqueryend":prevalign["queryend"], "alignquerystart":align["querystart"]})


