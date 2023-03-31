#!/usr/bin/env python3

import argparse, configparser
from argparse import RawTextHelpFormatter
import sys
from virsequel import Virsequel
from version import __version__

def main():
    description=f"""
	Identify viral sequences in rna seq data \n
    Version: {__version__}
	"""

    parser = argparse.ArgumentParser(description=description, formatter_class=RawTextHelpFormatter, add_help=False)

    parser.add_argument('-c', '--config', type=str, help='Path to configuration file')

    input_group = parser.add_argument_group('Required Arguments')
    input_group.add_argument('-1', '--read_1', type=str, required=False, help="Forward RNA fastq.gz data")
    input_group.add_argument('-2', '--read_2', type=str, required=False, help="Reverse RNA fastq.gz data")
    input_group.add_argument('-a', '--adapter', type=str, required=False, help="adapter file for bbduk trimming")
    input_group.add_argument('-o', '--output', type=str, required=False, help="output path")
    input_group.add_argument('-b', '--blastn', action='store_true', help='run blastn instead of diamond blastx')
    input_group.add_argument('-d', '--database', type=str, default='nr', help='the path of the database to use')
    
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
    optional.add_argument('-t', '--threads', type=str, required=False, default='8', help="Number of threads for non-spades processes")

    args=parser.parse_args()
    
    if args.config:
        config = configparser.ConfigParser()
        config.read(args.config)

        # Get options from the config file
        adapters = config.get('input', 'adapters')
        threads = config.get('input', 'threads')
        output = config.get('input', 'output')
        database = config.get('input', 'database')
        blastn=config.getboolean('input', 'blastn')
        sample_list = config['input']['sample_list'].split(',')
        
        # Check sample_list is provided in the config file
        if config.has_option('input', 'sample_list'):
            # Split the samples in the config file by "," and then strip any white space that might separate entries.
            sample_list = [s.strip() for s in config['input']['sample_list'].split(',')]
            for sample in sample_list:
                R1=sample
                R2=sample.replace('_R1', '_R2')
                pipeline = Virsequel(R1, R2, adapters, threads, output, database, blastn)
                pipeline.run_pipeline()
        else:
            print("Error: Sample list missing from config file")
            sys.exit(1)

    else:
        # Get options from the command line
        R1 = args.read_1
        R2 = args.read_2
        adapters = args.adapter
        threads = args.threads
        output = args.output
        database=args.database
        blastn=args.blastn

        if not R1 and R2:
            print("Error: You must provide forward and reverse reads")
            parser.print_help()
            sys.exit(1)
        elif not adapters:
            print("Error: An adapter file is needed for trimming")
            parser.print_help()
            sys.exit(1)
        elif not output:
            print("Error: An output path is required")
            parser.print_help()
            sys.exit(1)
        elif not database:
            print("Error: Please provide a path to a blast database")
            parser.print_help()
            sys.exit(1)
        else:
            pipeline = Virsequel(R1, R2, adapters,threads, output, database, blastn)
            pipeline.run_pipeline()

if __name__ == '__main__':
     main()