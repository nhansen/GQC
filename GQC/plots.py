import os
import glob
import re
import logging
import importlib.resources

logger = logging.getLogger(__name__)

def plot_benchmark_align_coverage(assemblyname:str, benchname:str, outputdir:str, benchparams:dict):
    rfile_res = importlib.resources.files("GQC").joinpath('BenchCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        genomefile = benchparams["genomeregions"]
        nlocfile = benchparams["nstretchregions"]
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + genomefile + " " + nlocfile
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_testassembly_align_coverage(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('TestCoveragePlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + outputdir + " " + resourcedir + " " + benchname
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_mononuc_accuracy(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('MononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_qv_score_concordance(assemblyname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('PlotAssemblyQualValueAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_svcluster_align_plots(assemblyname:str, benchname:str, outputdir:str, refobj, mode='bench', prefix=''):
    rfile_res = importlib.resources.files("GQC").joinpath('PlotChromAligns.R')

    chromalignbedfiles = glob.glob(outputdir + "/" + prefix + "*.clusters.bed")
    returnvalues = []
    with importlib.resources.as_file(rfile_res) as rfile:
        for chrombed in chromalignbedfiles:
            if mode == 'bench':
                chromosome = chrombed.replace(".clusters.bed", "")
                chromosome = re.sub(r".*/*clustered_aligns\.", "", chromosome)
                chromlength = refobj.get_reference_length(chromosome)
                plotcommand = "Rscript " + str(rfile) + " " + chrombed + " " + assemblyname + " " + benchname + " " + outputdir + " " + str(chromlength)
                logger.debug(plotcommand)
                returnvalues.append(os.system(plotcommand))
            elif mode == 'compare':
                chromosome = chrombed.replace(".clusters.bed", "")
                chromosome = re.sub(r".*/*clustered_aligns\.", "", chromosome)
                chromlength = refobj.get_reference_length(chromosome)
                plotcommand = "Rscript " + str(rfile) + " " + chrombed + " " + assemblyname + " " + benchname + " " + outputdir + " " + str(chromlength)
                logger.debug(plotcommand)
                returnvalues.append(os.system(plotcommand))
            else:
                logger.critical("Unknown mode passed to plot_svcluster_align_plots: " + str(mode))
                returnvalues.append(1)

    return returnvalues

def plot_sv_indel_profile_plot(assemblyname:str, benchname:str, outputdir:str, resourcedir:str, refobj):
    rfile_res = importlib.resources.files("GQC").joinpath('IndelProfile.R')

    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + benchname + " " + outputdir + " " 
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)

    return returnvalue

def plot_assembly_error_stats(assemblyname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('IndelLengthPlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + assemblyname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_error_stats(readsetname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('IndelLengthPlot.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_mononuc_stats(readsetname:str, genomename:str, outputdir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('ReadMononucAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + genomename + " " + outputdir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

def plot_read_qv_score_concordance(readsetname:str, benchname:str, outputdir:str, resourcedir:str):
    rfile_res = importlib.resources.files("GQC").joinpath('PlotAssemblyQualValueAccuracy.R')
    with importlib.resources.as_file(rfile_res) as rfile:
        plotcommand = "Rscript " + str(rfile) + " " + readsetname + " " + benchname + " " + outputdir + " " + resourcedir
        logger.info(plotcommand)
        returnvalue = os.system(plotcommand)
    return returnvalue

