import os
import re
import sys
import random
from tempfile import mkstemp
import pysam
from pysam import VariantFile
import pybedtools
import logging
from GQC import bedtoolslib

logger = logging.getLogger(__name__)

def create_vcf_and_hcbed_from_gvcf(gvcffile:str, args)->int:
    [vcffile, hcbedfile] = ['', '']

    return [vcffile, hcbedfile]

def create_constructed_genome(reffasta:str, vcffile:str, hcbedfile:str, outputdir:str, outputfiles:dict, args)->int:
    samplename = args.samplename

    print("About to construct haplotype specific VCF files")
    [hap1vcf, hap2vcf] = create_haplotype_specific_vcffiles(vcffile, outputdir, outputfiles, args)

    hap1mutfasta = outputfiles['hap1mutfasta']
    hap2mutfasta = outputfiles['hap2mutfasta']
    
    hap1chain = outputfiles['hap1chain']
    hap2chain = outputfiles['hap2chain']

    create_hap_specific_fastafile(reffasta, hap1vcf, samplename + ".hap1", hap1mutfasta, hap1chain, outputdir, outputfiles, args)
    create_hap_specific_fastafile(reffasta, hap2vcf, samplename + ".hap2", hap2mutfasta, hap2chain, outputdir, outputfiles, args)

    hap1hcbed = outputfiles['hap1hcbed']
    hap2hcbed = outputfiles['hap2hcbed']

    hap1unmapped = outputfiles['hap1unmapped']
    hap2unmapped = outputfiles['hap2unmapped']

    hap1lowqualbed = outputfiles['hap1lowqualbed']
    hap2lowqualbed = outputfiles['hap2lowqualbed']

    hap1maskedfasta = outputfiles['hap1maskedfasta']
    hap2maskedfasta = outputfiles['hap2maskedfasta']

    hap1chroms = args.hap1chroms.split(',')
    hap2chroms = args.hap2chroms.split(',')

    lift_and_apply_hq_beds(hcbedfile, hap1chain, hap2chroms, hap1hcbed, hap1unmapped, hap1lowqualbed, hap1mutfasta, hap1maskedfasta, outputdir, outputfiles, args)
    lift_and_apply_hq_beds(hcbedfile, hap2chain, hap1chroms, hap2hcbed, hap2unmapped, hap2lowqualbed, hap2mutfasta, hap2maskedfasta, outputdir, outputfiles, args)

    return 1

def construct_genome_from_reference(refobj, hapvcf:str, hapmutfasta:str, hapchain:str, outputdir:str, outputfiles:dict, args)->list:

    return [hapmutfasta, hapchain]

def create_haplotype_specific_vcffiles(diploidvcf:str, outputdir:str, outputfiles:dict, args)->list:
    samplename = args.samplename

    hap1vcf = outputfiles['hap1vcf']
    hap2vcf = outputfiles['hap2vcf']

    hap1chroms = args.hap1chroms.split(',')
    hap2chroms = args.hap2chroms.split(',')

    vcf_in = VariantFile(diploidvcf)  # auto-detect input format
    vcfheaderstring = str(vcf_in.header)
    matchexp = r"(#CHROM.*" + re.escape(samplename) + ")"
    subsexp1 = r"\1.hap1"
    subsexp2 = r"\1.hap2"
    hap1headerstring = re.sub(matchexp, subsexp1, vcfheaderstring)
    hap2headerstring = re.sub(matchexp, subsexp2, vcfheaderstring)

    hap1headerfd, hap1headerpath = mkstemp()
    with os.fdopen(hap1headerfd, 'w') as fh:
        fh.write(hap1headerstring)
    fh.close()
    hap1vcfheader_in = VariantFile(hap1headerpath)  # auto-detect input format
    os.remove(hap1headerpath)

    hap2headerfd, hap2headerpath = mkstemp()
    with os.fdopen(hap2headerfd, 'w') as fh:
        fh.write(hap2headerstring)
    fh.close()
    hap2vcfheader_in = VariantFile(hap2headerpath)  # auto-detect input format
    os.remove(hap2headerpath)

    vcf_hap1_out = VariantFile(hap1vcf, 'w', header=hap1vcfheader_in.header)
    vcf_hap2_out = VariantFile(hap2vcf, 'w', header=hap2vcfheader_in.header)

    for rec in vcf_in.fetch():
        recsamples = rec.samples.keys()
        hap1chrom = False
        hap2chrom = False
        if rec.contig not in hap2chroms:
            #logger.debug(rec.contig + " is not in hap2 chroms " + args.hap2chroms + " so including in hap 1")
            hap1chrom = True
        if rec.contig not in hap1chroms:
            #logger.debug(rec.contig + " is not in hap1 chroms " + args.hap1chroms + " so including in hap 2")
            hap2chrom = True

        randnum = random.random()
        if hap1chrom:
            allelenum = 1
            newrec = convert_geno_to_haploid(rec, vcf_hap1_out, samplename, samplename + ".hap1", allelenum, randnum)
            vcf_hap1_out.write(newrec)
        if hap2chrom:
            allelenum = 2
            newrec = convert_geno_to_haploid(rec, vcf_hap2_out, samplename, samplename + ".hap2", allelenum, randnum)
            vcf_hap2_out.write(newrec)

    vcf_hap1_out.close()
    vcf_hap2_out.close()

    pysam.tabix_index(hap1vcf, preset="vcf", force=True)
    pysam.tabix_index(hap2vcf, preset="vcf", force=True)

    return [hap1vcf, hap2vcf]

def convert_geno_to_haploid(vcfrecord, vcffile, oldsamplename, newsamplename, allelenum:int, randnum):
    # Create a new record:
    recsamples = vcfrecord.samples.keys()
    newrec = vcffile.new_record(contig=vcfrecord.contig, start=vcfrecord.start, stop=vcfrecord.stop, filter=vcfrecord.filter, alleles = vcfrecord.alleles)

    # is the original genotype here phased or unphased?
    oldsample = vcfrecord.samples[oldsamplename]
    if oldsample.phased:
        if allelenum == 1:
            newallele = oldsample['GT'][0]
        else:
            newallele = oldsample['GT'][-1]
        newrec.samples[newsamplename]['GT'] = (newallele, newallele)
        newrec.samples[newsamplename].phased = True
    else:
        if randnum <= 0.5:
            if allelenum == 1:
                newallele = oldsample['GT'][0]
            else:
                newallele = oldsample['GT'][-1]
        else:
            if allelenum == 1:
                newallele = oldsample['GT'][-1]
            else:
                newallele = oldsample['GT'][-0]
        newrec.samples[newsamplename]['GT'] = (newallele, newallele)
        newrec.samples[newsamplename].phased = False
    
    return newrec

def create_hap_specific_fastafile(reffasta:str, hapvcf:str, hapsamplename:str, hapmutfasta:str, hapchain:str, outputdir:str, outputfiles:dict, args)->int:
    
    bcftoolscmd = "bcftools consensus -c " + hapchain + " -f " + reffasta + " -s " + hapsamplename + " " + hapvcf + " > " + hapmutfasta
    returnvalue = os.system(bcftoolscmd)

    return returnvalue

def lift_and_apply_hq_beds(hcbedfile:str, hcchain:str, hapchroms:list, haphcbed:str, hapunmapped:str, haplcbed:str, hapfasta:str, maskedhapfasta:str, outputdir:str, outputfiles:dict, args)->list:

    # lift the HQ bed to the new haplotype with the bcftools chain file:
    hcname = re.sub(".*/", "", hcbedfile)
    hcout = re.sub(".bed", "", hcname)
    liftcmd = "liftOver " + hcbedfile + " " + hcchain + " " + haphcbed + " " + hapunmapped + " > " + outputdir + "/" + hcout + ".liftover.out"
    logger.debug(liftcmd)
    liftreturnvalue = os.system(liftcmd)

    # subtract lifted bed from lifted genome bed to find low-qual portion:
    hcintervals = pybedtools.BedTool(haphcbed)
    hccorrecthapintervals = []
    for hapinterval in hcintervals:
        if hapinterval.chrom not in hapchroms:
            hccorrecthapintervals.append(hapinterval)

    hapgenomeintervals = bedtoolslib.genomeintervals(hapfasta)
    lcintervals = bedtoolslib.subtractintervals(hapgenomeintervals, hccorrecthapintervals)
    lcintervals.saveas(haplcbed)

    fastamaskcmd = "bedtools maskfasta -fi " + hapfasta + " -bed " + haplcbed + " -fo " + maskedhapfasta
    returnvalue = os.system(fastamaskcmd)

    return returnvalue

