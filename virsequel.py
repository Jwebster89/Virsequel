#!/usr/bin/env python3

import os, sys, subprocess, logging, datetime

class Virsequel():
    def __init__(self, R1, R2, adapters, threads, output):
            # self.path=path
            self.R1=R1
            self.R2=R2
            self.adapters=adapters
            self.threads=threads
            self.output=output
            # self.sampleID=os.path.basename(os.path.splitext(self.R1)[0]).replace(".fastq", "").replace(".fq", "")

    def init_logs(self,sampleID):
        # create logs directory if it does not exist
        log_dir = os.path.join("logs")
        if not os.path.exists(log_dir):
            os.makedirs(log_dir)
        # create a log file with timestamp, sample ID, and step
        timestamp = datetime.datetime.now()
        log_file = os.path.join(log_dir, f"{timestamp.strftime('%Y-%m-%d_%H-%M')}.{sampleID}.log")
        # configure logging
        logging.basicConfig(filename=log_file, level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

         
    def trimming(self, R1, R2, adapters, threads, outdir):
        
        # check if input files end with the expected file extensions
        valid_extensions = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
        if not any(R1.endswith(ext) for ext in valid_extensions) or not any(R2.endswith(ext) for ext in valid_extensions):
            print("Error: input files must be in fastq format")
            logging.error("Error: input files must be in fastq format")
            sys.exit(1)

        # get sample ID from the filename of R1 and remove file extension
        sampleID=os.path.basename(os.path.splitext(R1)[0]).replace(".fastq", "").replace(".fq", "")
        if sampleID.endswith("_R1"):
            sampleID = sampleID[:-3]

        # initialize logs
        self.init_logs(sampleID)

        # create directory to store trimmed reads if it does not exist
        trim_dir=os.path.join(outdir,"trimmed_reads")
        if not os.path.exists(trim_dir):
            os.mkdir(trim_dir)

        
        # log the start of the trimming process
        print(f"Trimming sample {sampleID} with bbduk and saving files to {trim_dir}")
        logging.info(f"Trimming sample {sampleID} with bbduk and saving files to {trim_dir}")
        
        # create output filenames for trimmed reads and stats file
        output_R1 = os.path.join(outdir, "trimmed_reads", os.path.basename(R1).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        output_R2 = os.path.join(outdir, "trimmed_reads", os.path.basename(R2).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        stats_file = os.path.join(outdir, "trimmed_reads", sampleID.replace("_R1", "").replace(".fastq", "").replace(".fq", "") + ".stats.txt")

        # create bbduk command
        cmd = ['bbduk.sh', '-Xmx2g', f'in1={R1}', f'in2={R2}', f'out1={output_R1}', f'out2={output_R2}', f'stats={stats_file}', f'ref={adapters}', 'ktrim=r', 'k=23', 'mink=11', 'hdist=1', 'qtrim=rl', 'trimq=10', 'maq=10', 'minlen=50', 'threads=48', 'tpe', 'tbo']

        # log the bbduk command
        logging.info(f"Running command: {' '.join(cmd)}")

        if not os.path.isfile(output_R1):
            try:
                # run bbduk and capture output/error messages to logs
                subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
            except subprocess.CalledProcessError as error:
                # log errors and raise exception if subprocess fails
                logging.error(f"Error while running command: {' '.join(cmd)}")
                logging.error(f"Error message: {error.stderr}")
                raise
            print(f"Finished trimming sample {sampleID}")
            logging.info(f"Finished trimming sample {sampleID}")
        else:
            logging.info(f"trimmed data for {sampleID} already exists. Skipping.")


    def assembly(self, R1, R2, outdir):
        assembly_dir=os.path.join(outdir,"assemblies")
        if not os.path.exists(assembly_dir):
            os.mkdir(assembly_dir)
                
        sampleID=os.path.basename(os.path.splitext(R1)[0]).replace(".fastq", "").replace(".fq", "")
        if sampleID.endswith("_R1"):
            sampleID = sampleID[:-3]

        # log the start of the assembly process
        print(f"Assembling sample {sampleID} with metaspades and saving files to {assembly_dir}")
        logging.info(f"Assembling sample {sampleID} with metaspades and saving files to {assembly_dir}")

        # create input filenames and output location for spades
        input_1=os.path.join(outdir, "trimmed_reads", os.path.basename(R1).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        input_2=os.path.join(outdir, "trimmed_reads", os.path.basename(R2).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        spades_out=os.path.join(assembly_dir,sampleID+"_spades")

         # create metaspades command
        cmd= ["metaspades.py", "-1" , input_1, "-2", input_2,  "--only-assembler", "-t", "8", "-m", "58", "-k", "21,31,41,51,61,71,81,91,101", "-o", spades_out]

        # log the spades command
        logging.info(f"Running command: {' '.join(cmd)}")
        print(f"Running metaspades on {sampleID}")
        if not os.path.isfile(os.path.join(spades_out, "contigs.fasta")):
            try:
                # run metaspades and capture output/error messages to logs
                subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
            except subprocess.CalledProcessError as error:
                # log errors and raise exception if subprocess fails
                logging.error(f"Error while running command: {' '.join(cmd)}")
                logging.error(f"Error message: {error.stderr}")
                raise
            print(f"Finished assembling sample {sampleID}")
            logging.info(f"Finished arimming sample {sampleID}")
        else:
            logging.info(f"Assembly for {sampleID} already exists. Skipping.")

    def blast(self):
         pass
    
    def run_pipeline(self):
        self.trimming(self.R1, self.R2, self.adapters, self.threads, self.output)
        self.assembly(self.R1, self.R2, self.output)
