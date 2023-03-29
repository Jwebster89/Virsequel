#!/usr/bin/env python3

import argparse
from argparse import RawTextHelpFormatter
from virsequel import Virsequel
from version import __version__

def main():
    description=f"""
	Identify viral sequences in rna seq data \n
    Version: {__version__}
	"""

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-1', '--read_1', type=str, required=True, help="Forward RNA fastq.gz data")
    required.add_argument('-2', '--read_2', type=str, required=True, help="Reverse RNA fastq.gz data")
    required.add_argument('-a', '--adapter', type=str, required=True, help="adapter file for bbduk trimming")
    required.add_argument('-o', '--output', type=str, required=True, help="output path")
    
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
    optional.add_argument('-t', '--threads', type=str, required=False, default='8', help="Number of threads for non-spades processes")

    args=parser.parse_args()
    R1=args.read_1
    R2=args.read_2
    adapters=args.adapter
    threads=args.threads
    output=args.output
    
    pipeline = Virsequel(R1, R2, adapters,threads, output)
    pipeline.run_pipeline()

if __name__ == '__main__':
     main()