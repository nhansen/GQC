import shutil
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def create_output_directory(directory)->None:
    #directory = args.prefix
    path = Path(directory)
    logger.info("Creating directory " + directory + " for output")
    path.mkdir(exist_ok=True)

    return path.as_posix()

def name_output_files(args, outputdir:str)->dict:
    files = {}
    files["nonincludedbed"] = outputdir + "/nonincludedregions." + args.benchmark + ".bed"
    files["allexcludedbed"] = outputdir + "/excludedregions." + args.benchmark + ".bed"
    files["alignplotdir"] = outputdir + "/alignmentplots"
    files["alignplotprefix"] = outputdir + "/alignmentplots/" + args.assembly + ".clustered_aligns"
    files["testgenomebed"] = outputdir + "/genome." + args.assembly + ".bed"
    if not args.n_bedfile:
        files["testnbed"] = outputdir + "/nlocs." + args.assembly + ".bed"
    else:
        files["testnbed"] = None
    files["testnonnbed"] = outputdir + "/atgcseq." + args.assembly + ".bed"
    files["phasemarkerbed"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".hapmers.bed"
    files["mergedphasemarkerbed"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".hapmers.merged.bed"
    files["mergedmarkerbedwithscaffnames"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".scaffnames.hapmers.merged.bed"
    files["phaseblockbed"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".phasedscaffolds.bed"
    files["hmmphaseblockbed"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".hmmphasedscaffolds.bed"
    files["mergedhmmphaseblockbed"] = outputdir + "/" + args.assembly + "." + args.benchmark + ".merged.hmmphasedscaffolds.bed"
    files["trimmedphasedalignprefix"] = outputdir + "/" + args.assembly + "_vs_" + args.benchmark + ".trimmedphased"
    files["aligntomatbenchprefix"] = args.assembly + "_vs_" + args.benchmark + ".mat"
    files["aligntopatbenchprefix"] = args.assembly + "_vs_" + args.benchmark + ".pat"
    files["matalignedregions"] = outputdir + "/" + args.benchmark + ".mat.covered." + args.assembly + ".bed"
    files["patalignedregions"] = outputdir + "/" + args.benchmark + ".pat.covered." + args.assembly + ".bed"
    files["truthcovered"] = outputdir + "/" + args.assembly + ".benchcovered." + args.benchmark + ".bed"
    files["testmatcovered"] = outputdir + "/testmatcovered." + args.assembly + ".bed"
    files["testpatcovered"] = outputdir + "/testpatcovered." + args.assembly + ".bed"
    files["generalstatsfile"] = outputdir + "/" + args.assembly + ".generalstats.txt"
    files["contiglengths"] = outputdir + "/" + args.assembly + ".contiglengths.txt"
    files["scaffoldlengths"] = outputdir + "/" + args.assembly + ".scaffoldlengths.txt"
    files["mononucstatsfile"] = outputdir + "/" + args.assembly + ".mononucstats.txt"
    files["structdetailsfile"] = outputdir + "/" + args.assembly + ".structurestats.txt"
    files["structvariantbed"] = outputdir + "/" + args.assembly + ".svs.bed"
    files["structvariantsvcf"] = outputdir + "/" + args.assembly + ".svs.vcf"
    files["clusterlengths"] = outputdir + "/" + args.assembly + ".alignclusterlengths.txt"
    files["coveredmononucsfile"] = outputdir + "/" + args.assembly + ".coveredmononucs." + args.benchmark + ".bed"
    files["mononucswithvariantsfile"] = outputdir + "/" + args.assembly + ".mononucswithvariants." + args.benchmark + ".bed"
    files["bencherrortypebed"] = outputdir + "/" + args.assembly + ".errortype." + args.benchmark + ".bed"
    files["benchexcludederrortypebed"] = outputdir + "/" + args.assembly + ".excludederrors." + args.benchmark + ".bed"
    files["testerrortypebed"] = outputdir + "/errortype." + args.assembly + ".bed"
    files["coveredhetsitealleles"] = outputdir + "/" + args.benchmark + ".coveredhetalleles." + args.assembly + ".bed"
    files["snvstatsfile"] = outputdir + "/" + args.assembly + ".singlenucerrorstats.txt"
    files["indelstatsfile"] = outputdir + "/" + args.assembly + ".indelerrorstats.txt"
    files["qvstatsfile"] = outputdir + "/" + args.assembly + ".qvstats.txt"

    return files

def name_read_stats_files(args, outputdir:str)->dict:
    files = {}
    files["includedbedfile"] = outputdir + "/" + args.readsetname + ".includedregions." + args.benchmark + ".bed"
    files["mononucstatsfile"] = outputdir + "/" + args.readsetname + ".mononucstats.txt"
    files["mononuchistfile"] = outputdir + "/" + args.readsetname + ".mononuchist.txt"
    files["mononuccompositionfile"] = outputdir + "/" + args.readsetname + ".mononuccomposition.txt"
    files["includedmononucfile"] = outputdir + "/" + args.readsetname + ".includedmononucs.txt"
    files["mononucoverviewfile"] = outputdir + "/" + args.readsetname + ".mononucgeneralstats.txt"
    files["includedstrprefix"] = outputdir + "/" + args.readsetname + ".includedstrs"
    files["strstatsprefix"] = outputdir + "/" + args.readsetname + ".strstats"
    files["strlengthdiffaccuracyprefix"] = outputdir + "/" + args.readsetname + ".strlengthdiffaccuracy"
    files["strlengthaccuracyprefix"] = outputdir + "/" + args.readsetname + ".strlengthaccuracy"
    files["readerrorfile"] = outputdir + "/" + args.readsetname + ".readerrors.txt"
    files["errorstatsfile"] = outputdir + "/" + args.readsetname + ".generalstats.txt"
    files["snvstatsfile"] = outputdir + "/" + args.readsetname + ".singlenucerrorstats.txt"
    files["indelstatsfile"] = outputdir + "/" + args.readsetname + ".indelerrorstats.txt"
    files["qvstatsfile"] = outputdir + "/" + args.readsetname + ".errorqvstats.txt"

    return files

def name_constructgenome_output_files(args, outputdir:str)->dict:
    files = {}
    prefix = args.prefix
    files["hap1vcf"] = outputdir + "/" + prefix + ".hap1.vcf.gz"
    files["hap2vcf"] = outputdir + "/" + prefix + ".hap2.vcf.gz"
    files["hap1chain"] = outputdir + "/" + prefix + ".hap1.chain"
    files["hap2chain"] = outputdir + "/" + prefix + ".hap2.chain"
    files["hap1mutfasta"] = outputdir + "/" + prefix + ".hap1.mut.fasta"
    files["hap2mutfasta"] = outputdir + "/" + prefix + ".hap2.mut.fasta"
    files["hap1hcbed"] = outputdir + "/" + prefix + ".hap1.highconf.bed"
    files["hap2hcbed"] = outputdir + "/" + prefix + ".hap2.highconf.bed"
    files["hap1unmapped"] = outputdir + "/" + prefix + ".hap1.unmapped"
    files["hap2unmapped"] = outputdir + "/" + prefix + ".hap2.unmapped"
    files["hap1lowqualbed"] = outputdir + "/" + prefix + ".hap1.lowconf.bed"
    files["hap2lowqualbed"] = outputdir + "/" + prefix + ".hap2.lowconf.bed"
    files["hap1maskedfasta"] = outputdir + "/" + prefix + ".hap1.masked.fasta"
    files["hap2maskedfasta"] = outputdir + "/" + prefix + ".hap2.masked.fasta"

    return files

def name_compare_output_files(args, outputdir:str)->dict:
    files = {}
    #files[""] = outputdir + "/" + args.A + "_vs_" + args.B + ".bed"

    return files

