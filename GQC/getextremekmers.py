import sys
import os
import shutil
import re
import argparse
import logging
from GQC import kmers

logger = logging.getLogger(__name__)

def init_argparse() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        usage="%(prog)s [OPTION] [FILE]...",
        description="Find kmers with extreme ATGC composition in a FastK database"
    )
    parser.add_argument(
        "-v", "--version", action="version",
        version = f"{parser.prog} version 0.1.0"
    )
    parser.add_argument('--fastkdb', required=True, help='FastK database root')
    parser.add_argument('-p', '--prefix', type=str, required=False, default="hetsites", help='prefix to use in output filenames')
    parser.add_argument('--debug', action='store_true', required=False, help='print verbose output to log file for debugging purposes')

    return parser

def parse_arguments(args):
    parser = init_argparse()
    args = parser.parse_args(args)

    return args

def check_for_histex():
    if shutil.which("Tabex") is None:
        print("You don\'t seem to have Tabex in your path. Please install Tabex")
        logger.critical("You don\'t seem to have Tabex in your path. Please install Tabex")
        exit(1)
    return 0

def find_extreme_kmers(fastkdbroot:str, args):

    kmers.find_extreme_kmers(fastkdbroot)

    return 0

def main() -> None:

    args = parse_arguments(sys.argv[1:])

    logfile = args.prefix + ".log"
    logformat = '%(asctime)s %(message)s'
    if args.debug:
        logging.basicConfig(filename=logfile, level=logging.DEBUG, format=logformat)
        logger.info('Logging verbose output for debugging.')
    else:
        logging.basicConfig(filename=logfile, level=logging.INFO, format=logformat)

    check_for_histex()

    find_extreme_kmers(args.fastkdb, args)

if __name__ == "__main__":
    main()


