import re
import math
import logging
import pybedtools

logger = logging.getLogger(__name__)

def write_general_assembly_stats(refobj, queryobj, contigregions, gapregions, outputfiles, args)->dict:

    bmstats = {}

    generalstatspath = outputfiles["generalstatsfile"]
    scafflengthfile = outputfiles["scaffoldlengths"]
    contiglengthfile = outputfiles["contiglengths"]

    # total bases on benchmark haplotypes (for NG/LG calculation):
    hap1totalbases = 0
    hap2totalbases = 0
    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    for refcontig in refobj.references:
        hap1match = phap1.match(refcontig)
        hap2match = phap2.match(refcontig)
        if hap1match:
            hap1totalbases = hap1totalbases + refobj.get_reference_length(refcontig)
        if hap2match:
            hap2totalbases = hap2totalbases + refobj.get_reference_length(refcontig)

    bmstats['hap1totalbases'] = hap1totalbases # this is the total number of bases in the MATERNAL benchmark chroms
    bmstats['hap2totalbases'] = hap2totalbases # this is the total number of bases in the PATERNAL benchmark chroms
    avghaptotalbases = (hap1totalbases + hap2totalbases)/2.0
    bmstats['totalbases'] = hap1totalbases + hap2totalbases

    mincontiglength = args.mincontiglength
    with open(generalstatspath, "w") as gsfh:
        numscaffolds = queryobj.nreferences
        totalscaffoldbases = sum(queryobj.lengths)

        cumscaffbases = 0
        numscaffs = 0
        scaffold_n50 = 0
        scaffold_hap1_ng50 = 0
        scaffold_hap2_ng50 = 0
        scaffold_l50 = 0
        scaffold_hap1_lg50 = 0
        scaffold_hap2_lg50 = 0
        scaffold_hap1_aung = 0
        scaffold_hap2_aung = 0
        scaffold_lengths = queryobj.lengths

        scaffold_lengths.sort(reverse=True)
        with open(scafflengthfile, "w") as sfh:
            for scafflength in scaffold_lengths:
                cumscaffbases = cumscaffbases + scafflength
                numscaffs = numscaffs + 1
                if cumscaffbases > 0.5*totalscaffoldbases and scaffold_n50 == 0:
                    scaffold_n50 = scafflength
                    scaffold_l50 = numscaffs
                if cumscaffbases > 0.5*hap1totalbases and scaffold_hap1_ng50 == 0:
                    scaffold_hap1_ng50 = scafflength
                    scaffold_hap1_lg50 = numscaffs
                if cumscaffbases > 0.5*hap2totalbases and scaffold_hap2_ng50 == 0:
                    scaffold_hap2_ng50 = scafflength
                    scaffold_hap2_lg50 = numscaffs
                if hap1totalbases > 0:
                    scaffold_hap1_aung = scaffold_hap1_aung + scafflength*scafflength/hap1totalbases
                if hap2totalbases:
                    scaffold_hap2_aung = scaffold_hap2_aung + scafflength*scafflength/hap2totalbases
                perctotallength = int(1000.0*cumscaffbases/avghaptotalbases + 0.5)/10.0
                # At some point, may want to pull scaffold name into last field for use in annotating plots
                sfh.write(str(perctotallength) + "\t" + str(scafflength) + "\t" + str(cumscaffbases) + "\tNA\n")

        totalns = 0
        for gap in gapregions:
            gaplength = len(gap)
            totalns = totalns + gaplength

        numcontigs = len(contigregions)
        contig_n50 = 0
        contig_hap1_ng50 = 0
        contig_hap2_ng50 = 0
        contig_l50 = 0
        contig_hap1_lg50 = 0
        contig_hap2_lg50 = 0
        contig_hap1_aung = 0
        contig_hap2_aung = 0
        numlargecontigs = 0
        sizelist = []
        totalsize = 0
        for contig in contigregions:
            contiglength = len(contig)
            if contiglength >= mincontiglength:
                numlargecontigs = numlargecontigs + 1
                sizelist.append(contiglength)
                totalsize = totalsize + contiglength

        sizelist.sort(reverse=True)
        cumcontigbases = 0
        numcontigs = 0
        with open(contiglengthfile, "w") as cfh:
            for contigsize in sizelist:
                cumcontigbases = cumcontigbases + contigsize
                numcontigs = numcontigs + 1
                if cumcontigbases > 0.5*totalsize and contig_n50 == 0:
                    contig_n50 = contigsize
                    contig_l50 = numcontigs
                if cumcontigbases > 0.5*hap1totalbases and contig_hap1_ng50 == 0:
                    contig_hap1_ng50 = contigsize
                    contig_hap1_lg50 = numcontigs
                if cumcontigbases > 0.5*hap2totalbases and contig_hap2_ng50 == 0:
                    contig_hap2_ng50 = contigsize
                    contig_hap2_lg50 = numcontigs
                if hap1totalbases > 0:
                    contig_hap1_aung = contig_hap1_aung + contigsize*contigsize/hap1totalbases
                if hap2totalbases > 0:
                    contig_hap2_aung = contig_hap2_aung + contigsize*contigsize/hap2totalbases
                perctotallength = int(1000.0*cumcontigbases/avghaptotalbases + 0.5)/10.0
                cfh.write(str(perctotallength) + "\t" + str(contigsize) + "\t" + str(cumcontigbases) + "\tNA\n")
        
        bmstats['totalns'] = totalns

        bmstats['numscaffolds'] = numscaffolds
        bmstats['totalscaffoldbases'] = totalscaffoldbases
        bmstats['scaffoldn50'] = scaffold_n50
        bmstats['scaffoldhap1ng50'] = scaffold_hap1_ng50
        bmstats['scaffoldhap2ng50'] = scaffold_hap2_ng50
        bmstats['scaffoldl50'] = scaffold_l50
        bmstats['scaffoldhap1lg50'] = scaffold_hap1_lg50
        bmstats['scaffoldhap2lg50'] = scaffold_hap2_lg50
        bmstats['scaffoldhap1aung'] = scaffold_hap1_aung
        bmstats['scaffoldhap2aung'] = scaffold_hap2_aung

        bmstats['numcontigs'] = numcontigs
        bmstats['numlargecontigs'] = numcontigs
        bmstats['totallargecontigbases'] = totalsize
        bmstats['contign50'] = contig_n50
        bmstats['contighap1ng50'] = contig_hap1_ng50
        bmstats['contighap2ng50'] = contig_hap2_ng50
        bmstats['contigl50'] = contig_l50
        bmstats['contighap1lg50'] = contig_hap1_lg50
        bmstats['contighap2lg50'] = contig_hap2_lg50
        bmstats['contighap1aung'] = contig_hap1_aung
        bmstats['contighap2aung'] = contig_hap2_aung

        gsfh.write("Scaffolds:\n\n")
        gsfh.write("Number of scaffolds: " + str(numscaffolds) + "\n")
        gsfh.write("Total bases in scaffolds: " + str(totalscaffoldbases) + "\n")
        gsfh.write("Ns per kb in scaffolds: " + str(int(totalns*1000/totalscaffoldbases)) + "\n")
        gsfh.write("Scaffold N50: " + str(round(scaffold_n50/1000000, 3)) + "Mb\n")
        gsfh.write("Scaffold NG50: ")
        gsfh.write(str(round(scaffold_hap2_ng50/1000000, 3)) + "Mb (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("Scaffold L50: " + str(scaffold_l50) + "\n")
        gsfh.write("Scaffold LG50: ")
        gsfh.write(str(scaffold_hap2_lg50) + " (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("Scaffold auNG: ")
        gsfh.write(str(round(scaffold_hap2_aung/1000000, 3)) + "Mb (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("\nContigs:\n\n")
        gsfh.write("Number of contigs: " + str(numcontigs) + "\n")
        gsfh.write("Number of contigs >= " + str(mincontiglength) + ": " + str(numlargecontigs) + "\n")
        gsfh.write("Total bases in contigs >= " + str(mincontiglength) + ": " + str(totalsize) + "\n")
        gsfh.write("Contig N50: " + str(round(contig_n50/1000000, 3)) + "Mb\n")
        gsfh.write("Contig NG50: ")
        gsfh.write(str(round(contig_hap2_ng50/1000000, 3)) + "Mb (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("Contig L50: " + str(contig_l50) + "\n")
        gsfh.write("Contig LG50: ")
        gsfh.write(str(contig_hap2_lg50) + " (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("Contig auNG: ")
        gsfh.write(str(round(contig_hap2_aung/1000000, 3)) + "Mb (calculated with respect to the length of the paternal haplotype)\n")
        gsfh.write("\n")

    return bmstats

def write_merged_aligned_stats(refobj, queryobj, mergedtruthcoveredbed, mergedtestmatcoveredbed, mergedtestpatcoveredbed, bedfiles:dict, bmstats:dict, args)->dict:

    generalstatsfile = bedfiles["generalstatsfile"]

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    totalbenchcovered = 0

    totalrefaligned = 0
    numberrefaligns = 0

    hap1totalbases = bmstats['hap1totalbases'] # total MATERNAL bases in benchmark
    hap2totalbases = bmstats['hap2totalbases'] # total PATERNAL bases in benchmark

    if hap1totalbases > 0:
        hap1_nga50 = 0
        hap1_lga50 = 0
        hap1_nga90 = 0
        hap1_lga90 = 0
        hap1_aunga = 0
    if hap2totalbases > 0:
        hap2_nga50 = 0
        hap2_lga50 = 0
        hap2_nga90 = 0
        hap2_lga90 = 0
        hap2_aunga = 0

    for truthint in sorted(mergedtruthcoveredbed, key=lambda h: len(h), reverse=True):
        alignlength = len(truthint)
        totalrefaligned = totalrefaligned + alignlength
        numberrefaligns = numberrefaligns + 1

        if hap1totalbases > 0:
            if totalrefaligned >= 0.5*hap1totalbases and hap1_nga50 == 0:
                hap1_nga50 = alignlength
                hap1_lga50 = numberrefaligns
            if totalrefaligned >= 0.9*hap1totalbases and hap1_nga90 == 0:
                hap1_nga90 = alignlength
                hap1_lga90 = numberrefaligns
            hap1_aunga = hap1_aunga + alignlength*alignlength/hap1totalbases
        if hap2totalbases > 0:
            if totalrefaligned >= 0.5*hap2totalbases and hap2_nga50 == 0:
                hap2_nga50 = alignlength
                hap2_lga50 = numberrefaligns
            if totalrefaligned >= 0.9*hap2totalbases and hap2_nga90 == 0:
                hap2_nga90 = alignlength
                hap2_lga90 = numberrefaligns
            hap2_aunga = hap2_aunga + alignlength*alignlength/hap2totalbases

    matbenchcovered = 0
    patbenchcovered = 0
    for truthint in mergedtruthcoveredbed:
        [chrom, start, end, name] = truthint
        chrom = truthint.chrom
        start = int(truthint.start)
        end = int(truthint.stop)
        totalbenchcovered = totalbenchcovered + end - start
        if phap1.match(chrom):
            matbenchcovered = matbenchcovered + end - start
        if phap2.match(chrom):
            patbenchcovered = patbenchcovered + end - start

    longesttestalignment = 0
    totaltestmatcovered = 0
    totaltestpatcovered = 0
    for testint in mergedtestmatcoveredbed:
        alignlength = len(testint)
        totaltestmatcovered = totaltestmatcovered + alignlength
        if alignlength > longesttestalignment:
            longesttestalignment = alignlength

    for testint in mergedtestpatcoveredbed:
        alignlength = len(testint)
        totaltestpatcovered = totaltestpatcovered + alignlength
        if alignlength > longesttestalignment:
            longesttestalignment = alignlength

    bmstats["testmattotalcovered"] = totaltestmatcovered
    bmstats["testpattotalcovered"] = totaltestpatcovered
    bmstats["benchtotalcovered"] = totalbenchcovered
    if hap1totalbases > 0:
        bmstats["mataunga"] = hap1_aunga
    if hap2totalbases > 0:
        bmstats["pataunga"] = hap2_aunga

    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nAligned contig bases:\n\n")
        if hap1totalbases > 0:
            gsfh.write("NGA50 (for MATERNAL benchmark haplotype): " + str(round(hap1_nga50/1000000, 3)) + "Mb\n")
            gsfh.write("LGA50 (for MATERNAL benchmark haplotype): " + str(hap1_lga50) + "\n")
            gsfh.write("NGA90 (for MATERNAL benchmark haplotype): " + str(round(hap1_nga90/1000000, 3)) + "Mb\n")
            gsfh.write("LGA90 (for MATERNAL benchmark haplotype): " + str(hap1_lga90) + "\n")
            gsfh.write("auNGA (for MATERNAL benchmark haplotype): " + str(round(hap1_aunga/1000000, 3)) + "Mb\n")
        if hap2totalbases > 0:
            gsfh.write("NGA50 (for PATERNAL benchmark haplotype): " + str(round(hap2_nga50/1000000, 3)) + "Mb\n")
            gsfh.write("LGA50 (for PATERNAL benchmark haplotype): " + str(hap2_lga50) + "\n")
            gsfh.write("NGA90 (for PATERNAL benchmark haplotype): " + str(round(hap2_nga90/1000000, 3)) + "Mb\n")
            gsfh.write("LGA90 (for PATERNAL benchmark haplotype): " + str(hap2_lga90) + "\n")
            gsfh.write("auNGA (for PATERNAL benchmark haplotype): " + str(round(hap2_aunga/1000000, 3)) + "Mb\n")
        perctestmatcovered = int(totaltestmatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        perctestpatcovered = int(totaltestpatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        if bmstats['hap1totalbases'] > 0 or bmstats['hap2totalbases'] > 0:
            percbenchcovered = int(1000 * totalbenchcovered / (bmstats['hap1totalbases'] + bmstats['hap2totalbases']) + 0.5) / 10
        else:
            percbenchcovered = 'NA'
        if bmstats['hap1totalbases'] > 0:
            percmatbenchcovered = int(1000 * matbenchcovered / bmstats['hap1totalbases'] + 0.5) / 10
        else:
            percmatbenchcovered = 'NA'
        if bmstats['hap2totalbases'] > 0:
            percpatbenchcovered = int(1000 * patbenchcovered / bmstats['hap2totalbases'] + 0.5) / 10
        else:
            percpatbenchcovered = 'NA'
        gsfh.write("Total " + args.assembly + " bases in aligns to MAT chromosomes: " + str(totaltestmatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestmatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        gsfh.write("Total " + args.assembly + " bases in aligns to PAT chromosomes: " + str(totaltestpatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestpatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        gsfh.write("Longest " + args.assembly + " alignment length (in test assembly bases) to the benchmark: " + str(round(longesttestalignment/1000000, 3)) + "Mb\n")
        gsfh.write("Total " + args.benchmark + " covered: " + str(totalbenchcovered) + "/" + str(bmstats['totalbases']) + " (" + str(percbenchcovered) + "% of all MAT and PAT bases in benchmark)" + "\n" )
        gsfh.write("MAT " + args.benchmark + " covered: " + str(matbenchcovered) + "/" + str(bmstats['hap1totalbases']) + " (" + str(percmatbenchcovered) + "% of MAT bases in benchmark)" + "\n")
        gsfh.write("PAT " + args.benchmark + " covered: " + str(patbenchcovered) + "/" + str(bmstats['hap2totalbases']) + " (" + str(percpatbenchcovered) + "% of PAT bases in benchmark)" + "\n")

    return bmstats

def write_aligned_cluster_stats(outputfiles:dict, bmstats:dict, args)->int:
    structstatsfile = outputfiles['clusterlengths']

    alignclusters = bmstats['alignclusters']
    hap1totalbases = bmstats['hap1totalbases'] # total MATERNAL bases in benchmark
    hap2totalbases = bmstats['hap2totalbases'] # total PATERNAL bases in benchmark
    avghaptotalbases = (hap1totalbases + hap2totalbases)/2.0

    allclusters = []

    for refentry in alignclusters.keys():
        refclusters = alignclusters[refentry]
        for cluster in refclusters:
            clusteraligns = cluster["aligns"]
            clusterlow = None
            clusterhigh = None
            for align in clusteraligns:
                if clusterlow is None or align["targetstart"] < clusterlow:
                    clusterlow = align["targetstart"]
                if clusterhigh is None or align["targetend"] > clusterhigh:
                    clusterhigh = align["targetend"]

            logger.debug(refentry + " cluster " + str(clusterlow) + "-" + str(clusterhigh))
            cluster["spanlength"] = clusterhigh - clusterlow + 1
            cluster["refentry"] = refentry
            allclusters.append(cluster)

    allclusters.sort(key=lambda x:x["spanlength"], reverse=True)
    totallength = 0
    with open(structstatsfile, "w") as cfh:
        for cluster in allclusters:
            spanlength = cluster["spanlength"]
            totallength = totallength + spanlength
            perctotallength = int(1000.0*totallength/avghaptotalbases + 0.5)/10.0
            refentry = cluster["refentry"]
            logger.debug(str(spanlength) + "\t" + str(totallength) + "\t" + refentry)
            cfh.write(str(perctotallength) + "\t" + str(spanlength) + "\t" + str(totallength) + "\t" + refentry + "\n")

    return 0

def write_aligned_stats(refobj, queryobj, truthcoveredbed, bedfiles:dict, bmstats:dict, args)->dict:

    generalstatsfile = bedfiles["generalstatsfile"]

    phap1 = re.compile(r'.*MAT.*')
    phap2 = re.compile(r'.*PAT.*')
    totalbenchcovered = 0

    totalrefaligned = 0
    numberrefaligns = 0

    hap1totalbases = bmstats['hap1totalbases'] # total MATERNAL bases in benchmark

    hap1_nga50 = 0
    hap1_lga50 = 0
    hap1_nga90 = 0
    hap1_lga90 = 0
    hap1_aunga = 0
    for truthint in sorted(truthcoveredbed, key=lambda h: len(h), reverse=True):
        alignlength = len(truthint)
        totalrefaligned = totalrefaligned + alignlength
        numberrefaligns = numberrefaligns + 1

        if totalrefaligned >= 0.5*hap1totalbases and hap1_nga50 == 0:
            hap1_nga50 = alignlength
            hap1_lga50 = numberrefaligns
        if totalrefaligned >= 0.9*hap1totalbases and hap1_nga90 == 0:
            hap1_nga90 = alignlength
            hap1_lga90 = numberrefaligns
        hap1_aunga = hap1_aunga + alignlength*alignlength/hap1totalbases

    matbenchcovered = 0
    patbenchcovered = 0
    for truthint in truthcoveredbed:
        [chrom, start, end, name] = truthint
        chrom = truthint.chrom
        start = int(truthint.start)
        end = int(truthint.stop)
        totalbenchcovered = totalbenchcovered + end - start
        if phap1.match(chrom):
            matbenchcovered = matbenchcovered + end - start
        if phap2.match(chrom):
            patbenchcovered = patbenchcovered + end - start

    #longesttestalignment = 0
    #totaltestmatcovered = 0
    #totaltestpatcovered = 0
    #for testint in mergedtestmatcoveredbed:
        #alignlength = len(testint)
        #totaltestmatcovered = totaltestmatcovered + alignlength
        #if alignlength > longesttestalignment:
            #longesttestalignment = alignlength
#
    #for testint in mergedtestpatcoveredbed:
        #alignlength = len(testint)
        #totaltestpatcovered = totaltestpatcovered + alignlength
        #if alignlength > longesttestalignment:
            #longesttestalignment = alignlength
#
    #bmstats["testmattotalcovered"] = totaltestmatcovered
    #bmstats["testpattotalcovered"] = totaltestpatcovered
    bmstats["benchtotalcovered"] = totalbenchcovered
    bmstats["mataunga"] = hap1_aunga

    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nAligned contig bases:\n\n")
        gsfh.write("NGA50 (for MATERNAL benchmark haplotype): " + str(round(hap1_nga50/1000000, 3)) + "Mb\n")
        gsfh.write("LGA50 (for MATERNAL benchmark haplotype): " + str(hap1_lga50) + "\n")
        gsfh.write("NGA90 (for MATERNAL benchmark haplotype): " + str(round(hap1_nga90/1000000, 3)) + "Mb\n")
        gsfh.write("LGA90 (for MATERNAL benchmark haplotype): " + str(hap1_lga90) + "\n")
        gsfh.write("auNGA (for MATERNAL benchmark haplotype): " + str(round(hap1_aunga/1000000, 3)) + "Mb\n")
        #perctestmatcovered = int(totaltestmatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        #perctestpatcovered = int(totaltestpatcovered * 1000 / bmstats['totallargecontigbases'] + 0.5) / 10
        if bmstats['hap1totalbases'] > 0 or bmstats['hap2totalbases'] > 0:
            percbenchcovered = int(1000 * totalbenchcovered / (bmstats['hap1totalbases'] + bmstats['hap2totalbases']) + 0.5) / 10
        else:
            percbenchcovered = 'NA'
        if bmstats['hap1totalbases'] > 0:
            percmatbenchcovered = int(1000 * matbenchcovered / bmstats['hap1totalbases'] + 0.5) / 10
        else:
            percmatbenchcovered = 'NA'
        if bmstats['hap2totalbases'] > 0:
            percpatbenchcovered = int(1000 * patbenchcovered / bmstats['hap2totalbases'] + 0.5) / 10
        else:
            percpatbenchcovered = 'NA'
        #gsfh.write("Total " + args.assembly + " bases in aligns to MAT chromosomes: " + str(totaltestmatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestmatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        #gsfh.write("Total " + args.assembly + " bases in aligns to PAT chromosomes: " + str(totaltestpatcovered) + "/" + str(bmstats['totallargecontigbases']) + " (" + str(perctestpatcovered) + "% of total bases in contigs >= " + str(args.mincontiglength) + " bases)" + "\n")
        #gsfh.write("Longest " + args.assembly + " alignment length (in test assembly bases) to the benchmark: " + str(round(longesttestalignment/1000000, 3)) + "Mb\n")
        gsfh.write("Total " + args.benchmark + " covered: " + str(totalbenchcovered) + "/" + str(bmstats['totalbases']) + " (" + str(percbenchcovered) + "% of all MAT and PAT bases in benchmark)" + "\n" )
        gsfh.write("MAT " + args.benchmark + " covered: " + str(matbenchcovered) + "/" + str(bmstats['hap1totalbases']) + " (" + str(percmatbenchcovered) + "% of MAT bases in benchmark)" + "\n")
        gsfh.write("PAT " + args.benchmark + " covered: " + str(patbenchcovered) + "/" + str(bmstats['hap2totalbases']) + " (" + str(percpatbenchcovered) + "% of PAT bases in benchmark)" + "\n")

    return bmstats

def write_qv_stats(benchmark_stats:dict, alignedscorecounts:list, snverrorscorecounts:list, indelerrorscorecounts:list, bedfiles:dict, args):

    variantfile = bedfiles["bencherrortypebed"]
    totalerrors = 0
    totalindelerrors = 0
    totalsnverrors = 0
    totalphasingerrors = 0
    totalphasingindelerrors = 0
    totalphasingsnverrors = 0
    totalconsensuserrors = 0
    totalconsensusindelerrors = 0
    totalconsensussnverrors = 0

    with open(variantfile, "r") as vfh:
        errorline = vfh.readline()
        while errorline:
            errorline = errorline.rstrip()
            [chrom, start, end, name, score, strand, widestart, wideend, color, errortype, vartype, varname] = errorline.split("\t")
            namefields = name.split("_")
            refallele = namefields[-2]
            altallele = namefields[-1]
            refallelelength = len(refallele) # these lengths are wrong when allele is "*"--should replace
            altallelelength = len(altallele) # these lengths are wrong when allele is "*"--should replace

            totalerrors = totalerrors + 1
            #if refallele == "*" or altallele == "*" or refallelelength != altallelelength:
            if vartype == 'INDEL':
                totalindelerrors = totalindelerrors + 1
            else:
                totalsnverrors = totalsnverrors + 1
            if errortype == "PHASING":
                totalphasingerrors = totalphasingerrors + 1
                if refallele == "*" or altallele == "*" or refallelelength != altallelelength:
                    totalphasingindelerrors = totalphasingindelerrors + 1
                else:
                    totalphasingsnverrors = totalphasingsnverrors + 1
            elif errortype == "CONSENSUS":
                totalconsensuserrors = totalconsensuserrors + 1
                if vartype == 'INDEL':
                    totalconsensusindelerrors = totalconsensusindelerrors + 1
                else:
                    totalconsensussnverrors = totalconsensussnverrors + 1

            errorline = vfh.readline()

    totalassemblybasesinaligns = benchmark_stats["testmattotalcovered"] + benchmark_stats["testpattotalcovered"]
    if totalassemblybasesinaligns > 0:
        phaseerrorrate = totalphasingerrors / totalassemblybasesinaligns
        phasesnverrorrate = totalphasingsnverrors / totalassemblybasesinaligns
        phaseindelerrorrate = totalphasingindelerrors / totalassemblybasesinaligns
        consensuserrorrate = totalconsensuserrors / totalassemblybasesinaligns
        consensussnverrorrate = totalconsensussnverrors / totalassemblybasesinaligns
        consensusindelerrorrate = totalconsensusindelerrors / totalassemblybasesinaligns
    else:
        phaseerrorrate = 'NA'
        phasesnverrorrate = 'NA'
        phaseindelerrorrate = 'NA'
        consensuserrorrate = 'NA'
        consensussnverrorrate = 'NA'
        consensusindelerrorrate = 'NA'

    if consensuserrorrate != "NA":
        if consensuserrorrate > 0:
            consensusqv = int(-10 * math.log(consensuserrorrate, 10) + 0.5)
        else:
            consensusqv = 'Inf'

        if consensuserrorrate + phaseerrorrate > 0:
            qvwithphaseerrors = int(-10 * math.log(consensuserrorrate + phaseerrorrate, 10) + 0.5)
        else:
            qvwithphaseerrors = 'Inf'
    else:
        consensusqv = 'NA'
        qvwithphaseerrors = 'NA'

    snvtypestring = ""
    totaltypesnverrors = 0
    if totalassemblybasesinaligns > 0:
        for errortype in sorted(benchmark_stats["singlebasecounts"].keys()):
            typeerrors = benchmark_stats["singlebasecounts"][errortype]
            snvtypeerrorspermb = typeerrors/totalassemblybasesinaligns*1000000
            snvtypestring = snvtypestring + errortype + "\t" + str(typeerrors) + "\t" + str(snvtypeerrorspermb) + "\n"
            totaltypesnverrors = totaltypesnverrors + typeerrors
        snverrorspermb = int(totaltypesnverrors/totalassemblybasesinaligns*1000000)
    
        with open(bedfiles["snvstatsfile"], "w") as sfh:
            sfh.write(snvtypestring)
    else:
        logger.info("Not writing statistics regarding substitutions because no bases aligned!")

    if totalassemblybasesinaligns > 0:
        indellengthstring = ""
        totalhistindelerrors = 0
        for indellength in sorted(benchmark_stats["indellengthcounts"].keys()):
            indelcount = benchmark_stats["indellengthcounts"][indellength]
            indelerrorspermb = indelcount/totalassemblybasesinaligns*1000000
            indellengthstring = indellengthstring + str(indellength) + "\t" + str(indelcount) + "\t" + str(indelerrorspermb) + "\n"
            totalhistindelerrors = totalhistindelerrors + indelcount
        indelerrorspermb = int(totalhistindelerrors/totalassemblybasesinaligns*1000000)
    
        with open(bedfiles["indelstatsfile"], "w") as ifh:
            ifh.write(indellengthstring)
    else:
        logger.info("Not writing statistics regarding indels because no bases aligned!")

    generalstatsfile = bedfiles["generalstatsfile"]
    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nAssembly accuracy:\n\n")
        gsfh.write("Total errors in alignments of " + args.assembly + " to " + args.benchmark + ": " + str(totalerrors) + "\n")
        gsfh.write("Total substitution errors in alignments of " + args.assembly + " to " + args.benchmark + ": " + str(totalsnverrors) + "\n")
        gsfh.write("Total indel errors in alignments of " + args.assembly + " to " + args.benchmark + ": " + str(totalindelerrors) + "\n")
        gsfh.write("Errors where sequence matches the opposite haplotype: " + str(totalphasingerrors) + " (" + str(totalphasingsnverrors) + " subs. and " + str(totalphasingindelerrors) + " indels)" + "\n")
        gsfh.write("Errors where sequence doesn\'t match the opposite haplotype: " + str(totalconsensuserrors) + " (" + str(totalconsensussnverrors) + " subs. and " + str(totalconsensusindelerrors) + " indels)" + "\n")
        gsfh.write("Overall QV without wrong-haplotype errors: " + str(consensusqv) + "\n")
        gsfh.write("Overall QV including wrong-haplotype errors: " + str(qvwithphaseerrors) + "\n")

    if len(alignedscorecounts) > 0:
        with open(bedfiles["qvstatsfile"], "w") as qfh:
            for i in range(len(alignedscorecounts)):
                if i > 0 or alignedscorecounts[i] > 0:
                    alignedbases = alignedscorecounts[i]
                    snverrorbases = snverrorscorecounts[i]
                    indelerrorbases = indelerrorscorecounts[i]
                    qfh.write(str(i) + "\t" + str(snverrorbases) + "\t" + str(indelerrorbases) + "\t" + str(alignedbases) + "\n")


    generalstatsfile = bedfiles["generalstatsfile"]
    
    return 0

def write_mononuc_stats(mononucstats:dict, bedfiles:dict, benchmark_stats:dict, args):

    totalcoveredmononucs = len(mononucstats.keys())
    totalcorrect = 0
    totalwronglength = 0
    totalphaseswitches = 0
    totalcomplex = 0

    for variantname in mononucstats.keys():
        mononucdict = mononucstats[variantname]
        repeatedbase = mononucdict["base"]
        varianttype = mononucdict["type"]
        reflength = mononucdict["length"]
        newlength = mononucdict["assemblylength"]
        if newlength == -1:
            totalcomplex = totalcomplex + 1
        elif varianttype == "CORRECT":
            totalcorrect = totalcorrect + 1
        elif varianttype == "PHASING":
            totalphaseswitches = totalphaseswitches + 1
        elif varianttype == "CONSENSUS":
            totalwronglength = totalwronglength + 1

    generalstatsfile = bedfiles["generalstatsfile"]
    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nMononucleotide run accuracy:\n\n")
        gsfh.write("Total number of true assembly mononucleotide runs of ten or more bases covered by " + args.assembly + " alignments: " + str(totalcoveredmononucs) + "\n")
        if totalcoveredmononucs > 0:
            perccorrect = round(100*totalcorrect/totalcoveredmononucs, 3)
        else:
            perccorrect = 'NA'
        gsfh.write("Number of mononucleotide runs correct in assembly: " + str(totalcorrect) + " (" + str(perccorrect) + "%)" + "\n")
        if totalcoveredmononucs > 0:
            percphaseswitch = round(100*totalphaseswitches/totalcoveredmononucs, 3)
        else:
            percphaseswitch = 'NA'
        gsfh.write("Number of mononucleotide runs with assembly alleles matching the alternate haplotype: " + str(totalphaseswitches) + " (" + str(percphaseswitch) + "%)" + "\n")
        if totalcoveredmononucs > 0:
            percwronglength = round(100*totalwronglength/totalcoveredmononucs, 3)
        else:
            percwronglength = 'NA'
        gsfh.write("Number of mononucleotide runs with non-alternate alleles of fewer or more of the same base in the assembly: " + str(totalwronglength) + " (" + str(percwronglength) + "%)" + "\n")
        if totalcoveredmononucs > 0:
            perccomplex = round(100*totalcomplex/totalcoveredmononucs, 3)
        else:
            perccomplex = 'NA'
        gsfh.write("Number of mononucleotide runs with erroneous alleles other than extensions or contractions: " + str(totalcomplex) + " (" + str(perccomplex) + "%)" + "\n")

    return 0

def write_het_stats(bedfiles:dict, bmstats:dict, args):

    hetbedfile = bedfiles["coveredhetsitealleles"]
    num_switches = 0
    num_maternal = 0
    num_paternal = 0
    total_matching_hets = 0

    pmat = re.compile(r'.*MAT.*')
    ppat = re.compile(r'.*PAT.*')

    last_contig = ""
    last_hap = ""
    with open(hetbedfile, "r") as hfh:
        hetline = hfh.readline()
        while hetline:
            hetline = hetline.rstrip()
            [contig, start, end, name, allele, chrom, chromstart, chromend, widecontig,  widestart, wideend, phasetype] = hetline.split("\t")
            alignedhap = "NA"
            if pmat.match(chrom):
                alignedhap = "MAT"
            elif ppat.match(chrom):
                alignedhap = "PAT"
            else:
                logger.warning("Uncertain haplotype for benchmark chromosome " + chrom)
            curallele_hap = "NA"
            if (phasetype == "SAMEHAP" and alignedhap == "MAT") or (phasetype == "ALTHAP" and alignedhap == "PAT"):
                num_maternal = num_maternal + 1
                total_matching_hets = total_matching_hets + 1
                curallele_hap = "MAT"
            if (phasetype == "SAMEHAP" and alignedhap == "PAT") or (phasetype == "ALTHAP" and alignedhap == "MAT"):
                num_paternal = num_paternal + 1
                curallele_hap = "PAT"
                total_matching_hets = total_matching_hets + 1

            # record switch if we're on the same contig and the haplotype changed:
            if contig == last_contig and curallele_hap != last_hap:
                num_switches = num_switches + 1
            last_contig = contig
            last_hap = curallele_hap

            hetline = hfh.readline()


    bmstats["numhetswitches"] = num_switches
    bmstats["nummaternalhetalleles"] = num_maternal
    bmstats["numpaternalhetalleles"] = num_paternal
    bmstats["phaseswitchespermb"] = num_switches/bmstats["totallargecontigbases"]*1000000

    generalstatsfile = bedfiles["generalstatsfile"]
    with open(generalstatsfile, "a") as gsfh:
        gsfh.write("\nPhasing statistics (" + str(total_matching_hets) + " benchmark heterozygous sites with an aligned maternal or paternal assembly allele):\n\n")
        gsfh.write("Number of maternal alleles within alignments: " + str(num_maternal) + "\n")
        gsfh.write("Number of paternal alleles within alignments: " + str(num_paternal) + "\n")
        gsfh.write("Number of switches between maternal and paternal alleles within a contig alignment: " + str(num_switches) + "\n")
        gsfh.write("Switch rate per base within aligned sequence: " + str(bmstats["phaseswitchespermb"]) + " per Mb\n")

    return 0

def write_read_mononuc_stats(stats:dict, outputfiles:dict, args):

    histogramfile = outputfiles["mononuchistfile"]
    lengthcomposition = {}
    with open(histogramfile, "w") as msfh:
        msfh.write("BenchmarkLength\tReadLength\tNumCorrect\tNumHetAllele\tNumError\tNumComplex\n")
        for mononuclength in sorted(stats.keys()):
            lengthcomposition[mononuclength] = {'CORRECT':0, 'HET':0, 'INSERROR':0, 'DELERROR':0, 'COMPLEX':0, 'FLANKERROR':0}
            for readmononuclength in sorted(stats[mononuclength].keys()):
                numcorrect = 0
                if "CORRECT" in stats[mononuclength][readmononuclength].keys():
                    numcorrect = stats[mononuclength][readmononuclength]["CORRECT"]
                    lengthcomposition[mononuclength]['CORRECT'] = lengthcomposition[mononuclength]['CORRECT'] + numcorrect
                numhetallele = 0
                if "HET" in stats[mononuclength][readmononuclength].keys():
                    numhetallele = stats[mononuclength][readmononuclength]["HET"]
                    lengthcomposition[mononuclength]['HET'] = lengthcomposition[mononuclength]['HET'] + numhetallele
                numerror = 0
                if "LENGTHERROR" in stats[mononuclength][readmononuclength].keys():
                    numerror = stats[mononuclength][readmononuclength]["LENGTHERROR"]
                    if readmononuclength > mononuclength:
                        lengthcomposition[mononuclength]['INSERROR'] = lengthcomposition[mononuclength]['INSERROR'] + numerror
                    elif readmononuclength < mononuclength:
                        lengthcomposition[mononuclength]['DELERROR'] = lengthcomposition[mononuclength]['DELERROR'] + numerror
                numcomplexallele = 0
                if "COMPLEX" in stats[mononuclength][readmononuclength].keys():
                    numcomplexallele = stats[mononuclength][readmononuclength]["COMPLEX"]
                    lengthcomposition[mononuclength]['COMPLEX'] = lengthcomposition[mononuclength]['COMPLEX'] + numcomplexallele
                if "FLANKERROR" in stats[mononuclength][readmononuclength].keys():
                    numflankallele = stats[mononuclength][readmononuclength]["FLANKERROR"]
                    lengthcomposition[mononuclength]['FLANKERROR'] = lengthcomposition[mononuclength]['FLANKERROR'] + numflankallele
                msfh.write(str(mononuclength) + "\t" + str(readmononuclength) + "\t" + str(numcorrect) + "\t" + str(numhetallele) + "\t" + str(numerror) + "\t" + str(numcomplexallele) + "\n")

    mononuccompfile = outputfiles["mononuccompositionfile"]
    mononucoverviewfile = outputfiles["mononucoverviewfile"]
    totalcorrect = 0
    totalerror = 0
    totalcomplex = 0
    totalinsertions = 0
    totaldeletions = 0
    totalflankerror = 0
    with open(mononuccompfile, "w") as mcfh:
        mcfh.write("BenchmarkLength\tNumCorrect\tNumHetAllele\tNumInsError\tNumDelError\tNumComplex\n")
        for mononuclength in sorted(stats.keys()):
            mcfh.write(str(mononuclength) + "\t" + str(lengthcomposition[mononuclength]['CORRECT']) + "\t" + str(lengthcomposition[mononuclength]['HET']) + "\t" + str(lengthcomposition[mononuclength]['INSERROR']) + "\t" + str(lengthcomposition[mononuclength]['DELERROR']) + "\t" + str(lengthcomposition[mononuclength]['COMPLEX']) + "\n")
            totalcorrect = totalcorrect + lengthcomposition[mononuclength]['CORRECT']
            totalerror = totalerror + lengthcomposition[mononuclength]['INSERROR'] + lengthcomposition[mononuclength]['DELERROR']
            totalcomplex = totalcomplex + lengthcomposition[mononuclength]['COMPLEX']
            totalflankerror = totalflankerror + lengthcomposition[mononuclength]['FLANKERROR']
            totalinsertions = totalinsertions + lengthcomposition[mononuclength]['INSERROR']
            totaldeletions = totaldeletions + lengthcomposition[mononuclength]['DELERROR']

    with open(mononucoverviewfile, "w") as mofh:
        if totalcorrect + totalerror > 0:
            noncomplexaccuracy = str(int(totalcorrect/(totalcorrect + totalerror)*1000)/10) + "%"
        else:
            noncomplexaccuracy = "NA"
        if totalcorrect + totalerror + totalcomplex > 0:
            accuracywithcomplex = str(int(totalcorrect/(totalcorrect + totalerror + totalcomplex)*1000)/10) + "%"
        else:
            accuracywithcomplex = "NA"
        if totalcorrect + totalerror + totalflankerror > 0:
            accuracywithflankerror = str(int(totalcorrect/(totalcorrect + totalerror + totalflankerror)*1000)/10) + "%"
        else:
            accuracywithflankerror = "NA"
        if totalcorrect + totalerror + totalcomplex + totalflankerror > 0:
            accuracywithcomplexandflankerror = str(int(totalcorrect/(totalcorrect + totalerror + totalcomplex + totalflankerror)*1000)/10) + "%"
        else:
            accuracywithcomplexandflankerror = "NA"

        if totaldeletions > 0:
            overcallundercallratio = str(int(100*totalinsertions/totaldeletions)/100)
        else:
            overcallundercallratio = "NA"
        mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (not counting read alleles that match the alternate haplotype or are complex): " + noncomplexaccuracy + "\n")
        mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including complex): " + accuracywithcomplex + "\n")
        mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including flankerror): " + accuracywithflankerror + "\n")
        mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including complex and flankerror): " + accuracywithcomplexandflankerror + "\n")
        mofh.write("Ratio of overcalled to undercalled homopolymer runs: " + overcallundercallratio + "\n")
            
def write_read_str_stats(strlabel:str, stats:dict, outputfiles:dict, args):

    lengthaccuracyfile = outputfiles["strlengthaccuracyprefix"] + "." + strlabel + ".txt"
    lengthdiffaccfile = outputfiles["strlengthdiffaccuracyprefix"] + "." + strlabel + ".txt"

    correctbylength = {}
    complexbylength = {}
    flankbylength = {}
    lengtherrorbylength = {}
    hetbylength = {}
   
    # tally length diff reads
    lengthdiffreads = {}
    complexreads = 0
    flankerrorreads = 0
    hetreads = 0
    
    for refstrlength in sorted(stats.keys()):
        correctbylength[refstrlength] = 0
        lengtherrorbylength[refstrlength] = 0
        complexbylength[refstrlength] = 0
        flankbylength[refstrlength] = 0
        hetbylength[refstrlength] = 0
        for readstrlength in sorted(stats[refstrlength].keys()):
            if readstrlength != -1:
                lengthdiff = readstrlength - refstrlength
                if lengthdiff not in lengthdiffreads.keys():
                    lengthdiffreads[lengthdiff] = 0
            else:
                lengthdiff = -1
            if "CORRECT" in stats[refstrlength][readstrlength].keys():
                correctbylength[refstrlength] = correctbylength[refstrlength] + stats[refstrlength][readstrlength]['CORRECT']
                lengthdiffreads[lengthdiff] = lengthdiffreads[lengthdiff] + stats[refstrlength][readstrlength]['CORRECT']
            if "LENGTHERROR" in stats[refstrlength][readstrlength].keys():
                lengtherrorbylength[refstrlength] = lengtherrorbylength[refstrlength] + stats[refstrlength][readstrlength]['LENGTHERROR']
                lengthdiffreads[lengthdiff] = lengthdiffreads[lengthdiff] + stats[refstrlength][readstrlength]['LENGTHERROR']
            if "COMPLEX" in stats[refstrlength][readstrlength].keys():
                complexbylength[refstrlength] = complexbylength[refstrlength] + stats[refstrlength][readstrlength]['COMPLEX']
                complexreads = complexreads + stats[refstrlength][readstrlength]['COMPLEX']
            if "FLANKERROR" in stats[refstrlength][readstrlength].keys():
                flankbylength[refstrlength] = flankbylength[refstrlength] + stats[refstrlength][readstrlength]['FLANKERROR']
                flankerrorreads = flankerrorreads + stats[refstrlength][readstrlength]['FLANKERROR']
            if "HET" in stats[refstrlength][readstrlength].keys():
                hetbylength[refstrlength] = hetbylength[refstrlength] + stats[refstrlength][readstrlength]['HET']
                hetreads = hetreads + stats[refstrlength][readstrlength]['HET']

    with open(lengthaccuracyfile, "w") as lfh:
        lfh.write("RefLength\tCorrect\tLengthError\tComplexError\tFlankError\tHet\n")
        for refstrlength in sorted(correctbylength.keys()):
            lfh.write(str(refstrlength) + "\t" + str(correctbylength[refstrlength]) + "\t" + str(lengtherrorbylength[refstrlength]) + "\t" + str(complexbylength[refstrlength]) + "\t" + str(flankbylength[refstrlength]) + "\t" + str(hetbylength[refstrlength]) + "\n")

    with open(lengthdiffaccfile, "w") as dfh:
        dfh.write("LengthDiff\tType\tNumReads\n")
        for lengthdiff in sorted(lengthdiffreads.keys()):
            if lengthdiff == 0:
                dfh.write(str(lengthdiff) + "\tCORRECT\t" + str(lengthdiffreads[lengthdiff]) + "\n")
            else:
                dfh.write(str(lengthdiff) + "\tLENGTHERROR\t" + str(lengthdiffreads[lengthdiff]) + "\n")
        dfh.write("0\tCOMPLEX\t" + str(complexreads) + "\n")
        dfh.write("0\tFLANKERROR\t" + str(flankerrorreads) + "\n")
        dfh.write("0\tHET\t" + str(hetreads) + "\n")

    #lengthcomposition = {}
    #with open(histogramfile, "w") as sfh:
        #sfh.write("BenchmarkLength\tReadLength\tNumCorrect\tNumHetAllele\tNumError\tNumComplex\n")
        #for mononuclength in sorted(stats.keys()):
            #lengthcomposition[mononuclength] = {'CORRECT':0, 'HET':0, 'INSERROR':0, 'DELERROR':0, 'COMPLEX':0, 'FLANKERROR':0}
            #for readmononuclength in sorted(stats[mononuclength].keys()):
                #numcorrect = 0
                #if "CORRECT" in stats[mononuclength][readmononuclength].keys():
                    #numcorrect = stats[mononuclength][readmononuclength]["CORRECT"]
                    #lengthcomposition[mononuclength]['CORRECT'] = lengthcomposition[mononuclength]['CORRECT'] + numcorrect
                #numhetallele = 0
                #if "HET" in stats[mononuclength][readmononuclength].keys():
                    #numhetallele = stats[mononuclength][readmononuclength]["HET"]
                    #lengthcomposition[mononuclength]['HET'] = lengthcomposition[mononuclength]['HET'] + numhetallele
                #numerror = 0
                #if "LENGTHERROR" in stats[mononuclength][readmononuclength].keys():
                    #numerror = stats[mononuclength][readmononuclength]["LENGTHERROR"]
                    #if readmononuclength > mononuclength:
                        #lengthcomposition[mononuclength]['INSERROR'] = lengthcomposition[mononuclength]['INSERROR'] + numerror
                    #elif readmononuclength < mononuclength:
                        #lengthcomposition[mononuclength]['DELERROR'] = lengthcomposition[mononuclength]['DELERROR'] + numerror
                #numcomplexallele = 0
                #if "COMPLEX" in stats[mononuclength][readmononuclength].keys():
                    #numcomplexallele = stats[mononuclength][readmononuclength]["COMPLEX"]
                    #lengthcomposition[mononuclength]['COMPLEX'] = lengthcomposition[mononuclength]['COMPLEX'] + numcomplexallele
                #if "FLANKERROR" in stats[mononuclength][readmononuclength].keys():
                    #numflankallele = stats[mononuclength][readmononuclength]["FLANKERROR"]
                    #lengthcomposition[mononuclength]['FLANKERROR'] = lengthcomposition[mononuclength]['FLANKERROR'] + numflankallele
                #sfh.write(str(mononuclength) + "\t" + str(readmononuclength) + "\t" + str(numcorrect) + "\t" + str(numhetallele) + "\t" + str(numerror) + "\t" + str(numcomplexallele) + "\n")
#
    #mononuccompfile = outputfiles["mononuccompositionfile"]
    #mononucoverviewfile = outputfiles["mononucoverviewfile"]
    #totalcorrect = 0
    #totalerror = 0
    #totalcomplex = 0
    #totalinsertions = 0
    #totaldeletions = 0
    #totalflankerror = 0
    #with open(mononuccompfile, "w") as mcfh:
        #mcfh.write("BenchmarkLength\tNumCorrect\tNumHetAllele\tNumInsError\tNumDelError\tNumComplex\n")
        #for mononuclength in sorted(stats.keys()):
            #mcfh.write(str(mononuclength) + "\t" + str(lengthcomposition[mononuclength]['CORRECT']) + "\t" + str(lengthcomposition[mononuclength]['HET']) + "\t" + str(lengthcomposition[mononuclength]['INSERROR']) + "\t" + str(lengthcomposition[mononuclength]['DELERROR']) + "\t" + str(lengthcomposition[mononuclength]['COMPLEX']) + "\n")
            #totalcorrect = totalcorrect + lengthcomposition[mononuclength]['CORRECT']
            #totalerror = totalerror + lengthcomposition[mononuclength]['INSERROR'] + lengthcomposition[mononuclength]['DELERROR']
            #totalcomplex = totalcomplex + lengthcomposition[mononuclength]['COMPLEX']
            #totalflankerror = totalflankerror + lengthcomposition[mononuclength]['FLANKERROR']
            #totalinsertions = totalinsertions + lengthcomposition[mononuclength]['INSERROR']
            #totaldeletions = totaldeletions + lengthcomposition[mononuclength]['DELERROR']
#
    #with open(mononucoverviewfile, "w") as mofh:
        #if totalcorrect + totalerror > 0:
            #noncomplexaccuracy = str(int(totalcorrect/(totalcorrect + totalerror)*1000)/10) + "%"
        #else:
            #noncomplexaccuracy = "NA"
        #if totalcorrect + totalerror + totalcomplex > 0:
            #accuracywithcomplex = str(int(totalcorrect/(totalcorrect + totalerror + totalcomplex)*1000)/10) + "%"
        #else:
            #accuracywithcomplex = "NA"
        #if totalcorrect + totalerror + totalflankerror > 0:
            #accuracywithflankerror = str(int(totalcorrect/(totalcorrect + totalerror + totalflankerror)*1000)/10) + "%"
        #else:
            #accuracywithflankerror = "NA"
        #if totalcorrect + totalerror + totalcomplex + totalflankerror > 0:
            #accuracywithcomplexandflankerror = str(int(totalcorrect/(totalcorrect + totalerror + totalcomplex + totalflankerror)*1000)/10) + "%"
        #else:
            #accuracywithcomplexandflankerror = "NA"
#
        #if totaldeletions > 0:
            #overcallundercallratio = str(int(100*totalinsertions/totaldeletions)/100)
        #else:
            #overcallundercallratio = "NA"
        #mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (not counting read alleles that match the alternate haplotype or are complex): " + noncomplexaccuracy + "\n")
        #mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including complex): " + accuracywithcomplex + "\n")
        #mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including flankerror): " + accuracywithflankerror + "\n")
        #mofh.write("Percent accuracy of homopolymer runs of 10 or more bases (including complex and flankerror): " + accuracywithcomplexandflankerror + "\n")
        #mofh.write("Ratio of overcalled to undercalled homopolymer runs: " + overcallundercallratio + "\n")

def write_read_error_summary(stats:dict, outputfiles:dict):
    totalsnverrors = 0
    totalindelerrors = 0
    snvtypestring = ""
    indellengthstring = ""
    with open(outputfiles["errorstatsfile"], "w") as efh:
        efh.write("Total aligned bases: " + str(stats["totalalignedbases"]) + "\n")
        efh.write("Total clipped bases: " + str(stats["totalclippedbases"]) + "\n")
        errorspermb = int(stats["totalerrorsinaligns"]/stats["totalalignedbases"]*1000000)
        efh.write("Total errors in alignments: " + str(stats["totalerrorsinaligns"]) + " (" + str(errorspermb) + " error per megabase)" + "\n")
        for errortype in sorted(stats["singlebasecounts"].keys()):
            typeerrors = stats["singlebasecounts"][errortype]
            snvtypeerrorspermb = typeerrors/stats["totalalignedbases"]*1000000
            snvtypestring = snvtypestring + errortype + "\t" + str(typeerrors) + "\t" + str(snvtypeerrorspermb) + "\n"
            totalsnverrors = totalsnverrors + typeerrors
        snverrorspermb = int(totalsnverrors/stats["totalalignedbases"]*1000000)
        efh.write("Total substitution errors in alignments: " + str(totalsnverrors) + " (" + str(snverrorspermb) + " error per megabase)" + "\n")
        for indellength in sorted(stats["indellengthcounts"].keys()):
            indelcount = stats["indellengthcounts"][indellength]
            indelerrorspermb = indelcount/stats["totalalignedbases"]*1000000
            indellengthstring = indellengthstring + str(indellength) + "\t" + str(indelcount) + "\t" + str(indelerrorspermb) + "\n"
            totalindelerrors = totalindelerrors + indelcount
        indelerrorspermb = int(totalindelerrors/stats["totalalignedbases"]*1000000)
        efh.write("Total indel errors in alignments: " + str(totalindelerrors) + " (" + str(indelerrorspermb) + " error per megabase)" + "\n")

    with open(outputfiles["snvstatsfile"], "w") as sfh:
        sfh.write(snvtypestring)

    with open(outputfiles["indelstatsfile"], "w") as ifh:
        ifh.write(indellengthstring)

    if len(stats["alignedqualscorecounts"]) > 0:
        with open(outputfiles["qvstatsfile"], "w") as qfh:
            for i in range(len(stats["alignedqualscorecounts"])):
                alignedbases = stats["alignedqualscorecounts"][i]
                snverrorbases = stats["snverrorqualscorecounts"][i]
                indelerrorbases = stats["indelerrorqualscorecounts"][i]
                qfh.write(str(i) + "\t" + str(snverrorbases) + "\t" + str(indelerrorbases) + "\t" + str(alignedbases) + "\n")
    #if len(stats["positiontotalcounts"]) > 0:
        #with open(outputfiles["positionstatsfile"], "w") as pfh:
            ## i is the zero-based position in the read
            #for i in range(len(stats["positionsnvcounts"])):
                #alignedbases = stats["positiontotalcounts"][i]
                #indelerrorbases = stats["positionindelcounts"][i]
                #pfh.write(str(i+1) + "\t" + str(indelerrorbases) + "\tIndel\t" + str(alignedbases) + "\n")
                #for snvkey in stats["positionsnvcounts"][i].keys():
                    #snverrorbases = stats["positionsnvcounts"][i]
                    #pfh.write(str(i+1) + "\t" + str(snverrorbases) + "\t" + snvkey + "\t" + str(alignedbases) + "\n")

    return 0

def write_phase_block_stats(phaseblockints:pybedtools.bedtool.BedTool, outputfiles:dict, benchmark_stats:dict, args):

    return 0

def write_assembly_haplotype_stats(fastaobj, contigregions, gapregions, args)->dict:

    hapstats = {}

    return hapstats

