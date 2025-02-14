import pytest
import subprocess
import sys
import os
import pysam
from qualiffy import bench
from qualiffy import output
from qualiffy import seqparse
from qualiffy import alignparse
from qualiffy import mummermethods
from qualiffy import bedtoolslib

def test_configs():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'blah', '-r', 'blah', '-q', 'blah', '-p', 'blah'])
    configvals = bench.read_config_data(args)
    assert configvals['resourcedir'] == 'resources'

def test_checkforprogs():
    bench.check_for_bedtools()
    no_rscript = bench.check_for_R()

def test_createoutputdir():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'tests/test.sort.bam', '-r', 'tests/testbenchmark.fasta.gz', '-q', 'tests/testassembly.fasta.gz', '-p', 'tests/testrun'])
    configvals = bench.read_config_data(args)
    outputdir = output.create_output_directory(args.prefix)
    assert os.path.isdir(outputdir)

def test_writebedfiles():
    args = bench.parse_arguments(['-c', 'tests/testconfig.txt', '-b', 'tests/test.sort.bam', '-r', 'tests/testbenchmark.fasta.gz', '-q', 'tests/testassembly.fasta.gz', '-p', 'tests/testrun'])
    configvals = bench.read_config_data(args)

    refobj = pysam.FastaFile(args.reffasta)
    queryobj = pysam.FastaFile(args.queryfasta)
    outputfiles = output.name_output_files(args, args.prefix)
    seqparse.write_genome_bedfiles(queryobj, refobj, args, configvals, outputfiles, {})

    with open('tests/testrun/genome.test.bed', 'r') as gfh:
        fields = gfh.readline().rstrip().split()
        assert(int(fields[2]) == 9721)
    with open('tests/testrun/atgcseq.test.bed', 'r') as gfh:
        fields = gfh.readline().rstrip().split()
        assert(int(fields[2]) == 9721)

    alignobj = pysam.AlignmentFile(args.bam, "rb")
    aligndata = alignparse.read_bam_aligns(alignobj, args.minalignlength)
    rlis_aligndata = mummermethods.filter_aligns(aligndata, "target")

    assert(len(aligndata) == 1)

