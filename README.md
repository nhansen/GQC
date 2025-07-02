# GQC -- Genome Quality Checker

The GQC python package analyses a user-supplied alignment of a test assembly to a benchmark assembly (preferably from the same sample), and prints general statistics, BED-formatted regions regarding the alignments and discrepancies within them, along with PDF-formatted plots. In addition, it reports statistics regarding discrepancies between aligned sequencing read sets and a benchmark that can help to elucidate sequencing platform strengths and weaknesses. Example outputs of various assembly and read benchmarking using GQC with the [HG002v1.1](https://github.com/marbl/HG002) assembly are available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=T2T/HG002/assemblies/benchmarking/analyses/).

The program was written by Nancy Fisher Hansen, a staff scientist in the [Genome Informatics Section](https://genomeinformatics.github.io/), Center for Genomics and Data Science Research, National Human Genome Research Institute (NHGRI). Nancy can be reached at nhansen@mail.nih.gov.

## Table of contents
- [Install](#install)
- [Getting Started](#getting-started)
- [Program Outputs](#program-outputs)

# Install

## Dependencies

Running GQC's assembly benchmarking requires an installation of commit 38b07c2 or later of Gene Myers' [FASTK](https://github.com/thegenemyers/FASTK.git) package with its "KmerMap" and "FastK" commands in your path. For both assembly and read analyses, GQC uses R's Rscript command with [Bioconductor](https://www.bioconductor.org/) to create plots, and [bedtools](https://bedtools.readthedocs.io/en/latest/) to compare and merge intervals. If the "Rscript" command is not in a user's path, the program will complain, and perform all functions except plotting. If the "bedtools" command isn't in the user's path, the program will exit with an error. To use GQC's plotting functions, you will need to install the "stringr" package and the "karyoploteR" package, which are part of Bioconductor.

In addition, the program requires a set of files with data about the benchmark assembly you are comparing to. For the Q100 benchmark assembly hg002v1.1, a tarball of these files is available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz). Once downloaded, this tarball should be unpacked and the locations of its files should be included in the config file that is either in your current working directory or that you pass to the GQC program (see the section "Config file" for more details).

All other dependencies will be installed by the pip installer with the commands in the next section called "Local Installation". Feel free to post installation issues to the issues section of this github repository and we will attempt to address them promptly.

## Local Installation

Until GQC becomes available on PyPi and bioconda, the easiest way to use it is to install it locally. First clone this github repository:
```
git clone https://github.com/nhansen/GQC
cd GQC
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

In order to evaluate heterozygous sites, short tandem repeat run lengths, and other features of the genome benchmark, the GQC program needs specially formatted annotation files for the benchmark genomes. For hg002v1.1, these files are contained in a tarball available on [AWS](https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/polishing/HG002/v1.1/benchmark/resources/hg002v1.1.resources.tar.gz). The program reads the locations of these files from a config file, which GQC assumes, by default, is a file called "benchconfig.txt" in your working directory, but the location of this file can also be specified with the -c or --config option. If the config file is not accessible in one of these two ways, the program will complain and exit.

An example config file is located in the resource tarball and contains the necessary parameters and file names. Edit that config file to specify the full path for each of the resource files (the "resourcedir" should be the path to the entire directory), and you can use that file as your config file when running the GQC or readbench commands.

## Evaluating haploid and diploid assemblies

To evaluate assembly scaffolds or contigs, the "bench" command first maps the locations of haplotype-specific 40-mers from the benchmark with the assembly's FASTA file. It then uses your installed version of minimap2 to create and trim alignments of each phased assembly sequence block to the appropriate benchmark haplotype.

	bench -r <benchmark.fasta> -q <assembly.fasta> -p <prefix_for_output> -A <assembly_name> -B <benchmark_name>

For typical assemblies, the bench command will use about 64Gb of memory and around 4 hours run time on two processors. The command "GQC --help" will display information on other options available (e.g., to restrict regions of the genome examined, set minimum contig or alignment lengths for processing, etc.).

## Evaluating read sets

To report and plot statistics about discrepancies between a set of sequencing reads and a benchmark diploid genome, the program has a "readbench" command. First, the reads should be aligned to the diploid benchmark assembly with whatever aligner and parameters you feel are most accurate. The usage of the readbench command is

	readbench -b <reads_vs_benchmark.bam> -r <benchmark.fasta> -p <prefix_for_output> -B <benchmark_name> -R <readset_name>

Because it is evaluating more alignments than for an assembly evaluation, the readbench command takes longer to run than the GQC command. For this reason, it has a "--downsample" option which allows the user to pass a fraction between 0 and 1.0 that will cause read alignments to be randomly downsampled to include only that fraction of the alignments in its accuracy calculations. As with the GQC command, information about options can be obtained with "readbench --help".

## Comparing two assemblies

To compare two FASTA files for two different assemblies (a "query" assembly and a "reference" assembly) of the same genome (which need not be a benchmark genome), the program "assemblycompare" first phases the query assembly against the reference assembly using 40 basepair k-mers that are unique to one haplotype of the reference assembly. It then reports statistics for completeness of and discrepancies within alignments between query scaffolds and the appropriate haplotype of the reference assembly. The usage for "assemblycompare" is

	assemblycompare --q1fasta <queryhap1.fasta> --q2fasta <queryhap2.fasta> --r1fasta <refhap1.fasta> --r2fasta <refhap2.fasta> -p <prefix_for_output> -Q <query_assembly_name> -R <ref_assembly_name>

# Program Outputs

All GQC programs create an output directory named with the prefix passed to the program with the --prefix (or -p) option. Within that directory will be a general statistics file, BED files, and pdf-formatted plots.

## GQC assembly evaluation outputs

### General statistics file

A file called "<assembly_name>.generalstats.txt" will contain general statistics about the assembly. The fasta file is first split into contigs anywhere the sequence contains at least 10 consecutive Ns (this value can be modified with the option --minns), and the number of scaffolds and contigs, the total bases within them, and statistics like N50/L50, NG50/LG50, and auNG are reported.

Then, the program evaluates the alignments of the phased scaffold segments to their corresponding benchmark genome haplotype, and reports similar statistics (NGA50/LGA50/auNGA) for these alignments, as well as the total number of bases aligned to each of the two benchmark haplotypes.

For each of the benchmark heterozygote sites provided in the benchmark's heterozygote bed file, the program will determine which of the two parental benchmark alleles is present in any of the haplotype alignments that cover the site. Then, for individual contig alignments, it will tally the number of times heterozygous sites within the alignment switch from maternal to paternal or vice-versa. From these switches and the lengths of the alignments, a switch rate is calculated and reported.

In evaluating accuracy, GQC tallies the number of discrepancies within primary alignments, and determines whether each represents the alternate allele of a heterozygous site. It then reports numbers of substitution and indel discrepancies and the number of these that match or don't match the alternate haplotype. From the total number of discrepancies and the total number of aligned bases, it reports a phred-scaled quality value (QV).

For each mono-, di-, tri-, and tetra- homopolymer run in the benchmark's short tandem repeat BED files, the bench program examines the assembly alignments that intersect the STR region. If the alignment has no discrepancies, it is considered "correct". If not and the assembly has a same-base run of a different length, its length is tallied, and if it has a different sequence, it is tallied as "complex". In the general statistics file, these categories are reported as "correct", "fewer or more of the same base", or "erroneous alleles other than extensions or contractions".

### BED files

The BED-formatted files produced by GQC include the following:

* excludedregions.<benchmark_name>.bed - regions that were excluded from analysis, either because they were in the config file's excluded regions, were not in regions specified with --includefile, and or were excluded with the --excludefile option
* nlocs.<assemblyhaplotype_name>.bed - regions with 10 or more Ns in the assembly haplotype (10 threshold can be modified with the --minns option)
* <assemblyhaplotype_name>.benchcovered.<benchmark_name>.bed - regions in the benchmark assembly covered alignments of assembly haplotype scaffolds
* testmatcovered.<assemlyhaplotype_name>.bed/testpatcovered.<assemblyhaplotype>.bed - regions in the assembly haplotype that align to maternal/paternal chromosomes in the benchmark
* <assemblyhaplotype_name>.errortype.<benchmark_name>.bed - benchmark locations of discrepancies in alignments, along with their locations in the test assembly haplotype
* <assemblyhaplotype_name>.mononucswithvariants.<benchmark_name>.bed - locations of benchmark mononucleotides covered by assembly alignments, with discrepancies

### Plots

Plot title names for the test assembly haplotype and the benchmark assembly are the names passed with the -A and -B options, respectively. The PDF-formatted plots currently produced by GQC are the following:

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
