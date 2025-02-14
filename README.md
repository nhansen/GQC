# qualiffy

The qualiffy python package analyses a user-supplied alignment of a test assembly to a benchmark assembly (preferably from the same sample), and prints general statistics, BED-formatted regions regarding the alignments and discrepancies within them, along with PDF-formatted plots.

The program was written by Nancy Fisher Hansen, a staff scientist in the Genome Informatics Section at the National Human Genome Research Institute (NHGRI). Nancy can be reached at nhansen@mail.nih.gov.

## Table of contents
- [Install](#install)
- [Getting Started](#getting-started)
- [Program Outputs](#program-outputs)

# Install

## Dependencies

This program uses R's Rscript command with [Bioconductor](https://www.bioconductor.org/) to create plots, and [bedtools](https://bedtools.readthedocs.io/en/latest/) to compare and merge intervals. If the "Rscript" command is not in a user's path, the program will complain, and perform all functions except plotting. If the "bedtools" command isn't in the user's path, the program will exit with an error. To use qualiffy's plotting functions, you will need to install the "stringr" package and the "karyoploteR" package, which is part of Bioconductor.

In addition, the program requires various files with data about the benchmark assembly you are comparing to. For the Q100 benchmark assembly hg002v1.1, a tarball of these files is available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz). Once downloaded, the tarball should be unpacked and the locations of its files should be included in the config file passed to the program (see the section "Config file" for more details).

All other dependencies will be installed by the pip installer with the commands in the next section called "Local Installation". Feel free to post installation issues to the issues section of this github repository.

## Local Installation

Until qualiffy is available on PyPi and bioconda, the easiest way to use it is to install it locally. First clone this github repository:
```
git clone https://github.com/nhansen/qualiffy
cd qualiffy
```

Create a virtual environment for the project:
```
python3 -m venv venv
source venv/bin/activate
```

Finally use python's pip installer to install and test a development copy for yourself to run:
```
python3 -m pip install -e .
pytest
```
# Getting Started

## Config file

In order to evaluate heterozygous sites, mononucleotide run lengths, and other features of the benchmark, the qualiffy program needs specially formatted annotation files for the benchmark genomes. For hg002v1.1, these files are contained in a tarball available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz). The program reads the locations of these files from a config file, which qualiffy assumes, by default, is a file called "benchconfig.txt" in your working directory. The location of this file can also be specified with the -c or --config option. If the config file is not accessible in one of these two ways, the program will complain and exit.

An example config file is located in the resource tarball and contains the necessary parameters and file names. Edit that config file to specify the full path for each of the resource files (the "resourcedir" should be the path to the entire directory), and you can use that file as your config file when running the qualiffy or readbench commands.

## Evaluating haplotypes of diploid assemblies

At this point, the qualiffy program is written to evaluate a BAM-formatted file of alignments of a single haplotype of a diploid assembly to a diploid benchmark genome (e.g., hg002v1.1). We recommend aligning each haplotype of the diploid assembly separately to the diploid benchmark using minimap2 with the "-x asm5" preset. To evaluate each haplotype the program usage is

	qualiffy assemblybench -b <assemblyhaplotype_vs_benchmark.bam> -r <benchmark.fasta> -q <assemblyhaplotype.fasta> -p <prefix_for_output> -A <assemblyhaplotype_name> -B <benchmark_name>

For typical assemblies, the assemblybench will use about 16Gb of memory and around 15 minutes of CPU time. The command "qualiffy --help" will display information on other options available (e.g., to restrict regions of the genome examined, set minimum contig or alignment lengths for processing, etc.).

## Evaluating read sets

To report and plot statistics about discrepancies between a set of sequencing reads and a benchmark diploid genome, the program has a "readbench" command. First, the reads should be aligned to the diploid benchmark assembly with whatever aligner and parameters you feel are most accurate. The usage of the readbench command is

	qualiffy readbench -b <reads_vs_benchmark.bam> -r <benchmark.fasta> -p <prefix_for_output> -B <benchmark_name> -R <readset_name>

Because it is evaluating more alignments than for an assembly evaluation, the readbench command takes longer to run than the qualiffy command. For this reason, it has a "--downsample" option which allows the user to pass a fraction between 0 and 1.0 that will cause read alignments to be randomly downsampled to include only that fraction of the alignments in its accuracy calculations. As with the qualiffy command, information about options can be obtained with "readbench --help".

# Program Outputs

All qualiffy programs create an output directory named with the prefix passed to the program with the --prefix (or -p) option. Within that directory will be a general statistics file, BED files, and pdf-formatted plots.

## qualiffy assembly evaluation outputs

### General statistics file

A file called "<assemblyhaplotype_name>.generalstats.txt" will contain general statistics about the assembly haplotype. The haplotype fasta file is first split into contigs anywhere the sequence contains at least 10 consecutive Ns (this value can be modified with the option --minns), and the number of scaffolds and contigs, the total bases within them, and statistics like N50/L50, NG50/LG50, and auNG are reported.

Then, the program evaluates the user-supplied alignments of the haplotype to the benchmark genome, and reports similar statistics (NGA50/LGA50/auNGA) for these alignments, as well as the total number of haplotype bases aligned to each of the two benchmark haplotypes.

For each of the benchmark heterozygote sites provided in the benchmark's heterozygote bed file, the program will determine which of the two parental benchmark alleles is present in any of the haplotype alignments that cover the site. Then, for individual contig alignments, it will tally the number of times heterozygous sites within the alignment switch from maternal to paternal or vice-versa. From these switches and the lengths of the alignments, a switch rate is calculated and reported.

In evaluating accuracy, qualiffy tallies the number of discrepancies within primary alignments, and determines whether each represents the alternate allele of a heterozygous site. It then reports numbers of substitution and indel discrepancies and the number of these that match or don't match the alternate haplotype. From the total number of discrepancies and the total number of aligned bases, it reports a phred-scaled quality value (QV).

For each homopolymer run in the benchmark's mononucleotide run file, the program examines the assembly alignments that intersect the homopolymer region. If the alignment has no discrepancies, it is considered "correct". If not and the assembly has a same-base run of a different length, its length is tallied, and if it has a different sequence, it is tallied as "complex". In the general statistics file, these categories are reported as "correct", "fewer or more of the same base", or "erroneous alleles other than extensions or contractions".

### BED files

The BED-formatted files produced by qualiffy include the following:

* excludedregions.<benchmark_name>.bed - regions that were excluded from analysis, either because they were in the config file's excluded regions, were not in regions specified with --includefile, and or were excluded with the --excludefile option
* nlocs.<assemblyhaplotype_name>.bed - regions with 10 or more Ns in the assembly haplotype (10 threshold can be modified with the --minns option)
* <assemblyhaplotype_name>.benchcovered.<benchmark_name>.bed - regions in the benchmark assembly covered alignments of assembly haplotype scaffolds
* testmatcovered.<assemlyhaplotype_name>.bed/testpatcovered.<assemblyhaplotype>.bed - regions in the assembly haplotype that align to maternal/paternal chromosomes in the benchmark
* <assemblyhaplotype_name>.errortype.<benchmark_name>.bed - benchmark locations of discrepancies in alignments, along with their locations in the test assembly haplotype
* <assemblyhaplotype_name>.mononucswithvariants.<benchmark_name>.bed - locations of benchmark mononucleotides covered by assembly alignments, with discrepancies

### Plots

Plot title names for the test assembly haplotype and the benchmark assembly are the names passed with the -A and -B options, respectively. The PDF-formatted plots currently produced by qualiffy are the following:

* <assemblyhaplotype_name>.benchcovered.<benchmark_name>.pdf - a karyotype plot of the benchmark diploid genome with maternal chromosomes on top and paternal chromosomes on bottom. Chromosomes are colored in wherever they are covered by an alignment of the tested assembly haplotype
* <assemblyhaplotype_name>.benchcoveredwitherrors.<benchmark_name>.pdf - same as the "benchcovered" plot, but with a wiggle plot of locations of discrepancies within alignments
* <assemblyhaplotype_name>.testcovered.<benchmark_name>.pdf - karyotype plot of test assembly scaffolds, colored to show where they align to maternal and paternal benchmark chromosomes
* <assemblyhaplotype_name>.indelsizestats.pdf - histogram of sizes of insertions and deletions within alignments of the test assembly haplotype to the benchmark genome
* <assemblyhaplotype_name>.mononuc_accuracy.<benchmark_name>.pdf - percent of mononucleotide runs accurate in the test assembly haplotype, plotted by mononucleotide length

## readbench sequencing read evaluation outputs

### General statistics file

A file called "<readset_name>.generalstats.txt" will report the total number of aligned read bases, the total number of clipped read bases, the total number of discrepancies within alignments (with a rate of discprepancies per aligned megabase), and the breakdown of these discrepancies into substitution and indel changes.

### Tab-delimited statistics files

* <readset_name>.singlenucerrorstats.txt - strand-specific nucleotide changes with the number observed and the rate they occur per aligned megabase. The first base in the reported change is the benchmark base (complemented if the read aligns along the reverse strand) and the second is the base reported within the read
* <readset_name>.indelerrorstats.txt - a tab-delimited file with the size of the observed insertion or deletion (negative is when bases from the benchmark are deleted in the read, positive when there are inserted bases in the read), the number observed, and the number observed per aligned megabase
* <readset_name>.readerrors.txt - a BED-formatted file with the locations of all tallied discrepancies between reads and the benchmark genome (keep in mind that when the --downsample option is used, these will not include all errors in all reads, just the ones in alignments that pass the downsampling threshold)

### Plots

* <readset_name>.indelsizestats.pdf - a histogram of the observed rates of indel discrepancies per aligned megabase of read
