#!/usr/bin/env python3

import os, sys, subprocess, logging, datetime

class Virsequel():
    def __init__(self, R1, R2, adapters, threads, output, db, blastn):
            # self.path=path
            self.R1=R1
            self.R2=R2
            self.adapters=adapters
            self.threads=threads
            self.output=output
            self.db=db
            self.blastn=blastn
            self.sampleID=os.path.basename(os.path.splitext(self.R1)[0]).replace(".fastq", "").replace(".fq", "")
            if self.sampleID.endswith("_R1"):
                self.sampleID = self.sampleID[:-3]
            
            # create logs directory if it does not exist
            log_dir = os.path.join("logs")
            if not os.path.exists(log_dir):
                os.makedirs(log_dir)
            # create a log file with timestamp, sample ID, and step
            timestamp = datetime.datetime.now()
            log_file = os.path.join(log_dir, f"{timestamp.strftime('%Y-%m-%d_%H-%M')}.{self.sampleID}.log")

            # configure logging
            self.logger = logging.getLogger(self.sampleID)
            self.logger.setLevel(logging.INFO)
            formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)

        
    def trimming(self, R1, R2, adapters, threads, outdir):
        print("STEP 1: Trimming")
        # logger = logging.getLogger(self.sampleID)
        # check if input files end with the expected file extensions
        valid_extensions = [".fq", ".fastq", ".fq.gz", ".fastq.gz"]
        if not any(R1.endswith(ext) for ext in valid_extensions) or not any(R2.endswith(ext) for ext in valid_extensions):
            print("Error: input files must be in fastq format")
            self.logger.error("Error: input files must be in fastq format")
            sys.exit(1)

        # create directory to store trimmed reads if it does not exist
        trim_dir=os.path.join(outdir,"trimmed_reads")
        if not os.path.exists(trim_dir):
            os.mkdir(trim_dir)

        # create output filenames for trimmed reads and stats file
        output_R1 = os.path.join(outdir, "trimmed_reads", os.path.basename(os.path.splitext(R1)[0]).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        output_R2 = os.path.join(outdir, "trimmed_reads", os.path.basename(os.path.splitext(R2)[0]).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        stats_file = os.path.join(outdir, "trimmed_reads", self.sampleID.replace("_R1", "").replace(".fastq", "").replace(".fq", "") + ".stats.txt")

        # create bbduk command
        cmd = ['bbduk.sh', '-Xmx2g', f'in1={R1}', f'in2={R2}', f'out1={output_R1}', f'out2={output_R2}', f'stats={stats_file}', f'ref={adapters}', 'ktrim=r', 'k=23', 'mink=11', 'hdist=1', 'qtrim=rl', 'trimq=10', 'maq=10', 'minlen=50', 'threads=48', 'tpe', 'tbo']
        
        # log the start of the trimming process
        print(f"Trimming sample {self.sampleID} with bbduk and saving files to {trim_dir}")
        self.logger.info(f"Trimming sample {self.sampleID} with bbduk and saving files to {output_R1} and {output_R2}")
        # log the bbduk command
        self.logger.info(f"Running command: {' '.join(cmd)}")

        if not os.path.isfile(output_R1):
            try:
                # run bbduk and capture output/error messages to logs
                subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
            except subprocess.CalledProcessError as error:
                # log errors and raise exception if subprocess fails
                self.logger.error(f"Error while running command: {' '.join(cmd)}")
                self.logger.error(f"Error message: {error.stderr}")
                raise
            print(f"Finished trimming sample {self.sampleID}")
            self.logger.info(f"Finished trimming sample {self.sampleID}")
        else:
            self.logger.info(f"trimmed data for {self.sampleID} already exists. Skipping.")


    def assembly(self, R1, R2, outdir):
        print("\nSTEP 2: Assembly")
        # logger = logging.getLogger(self.sampleID)
        assembly_dir=os.path.join(outdir,"assemblies")
        if not os.path.exists(assembly_dir):
            os.mkdir(assembly_dir)
                
        # sampleID=os.path.basename(os.path.splitext(R1)[0]).replace(".fastq", "").replace(".fq", "")
        if self.sampleID.endswith("_R1"):
            self.sampleID = self.sampleID[:-3]

        # log the start of the assembly process
        print(f"Assembling sample {self.sampleID} with metaspades and saving files to {assembly_dir}")
        self.logger.info(f"Assembling sample {self.sampleID} with metaspades and saving files to {assembly_dir}")

        # create input filenames and output location for spades
        input_1=os.path.join(outdir, "trimmed_reads", os.path.basename(os.path.splitext(R1)[0]).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        input_2=os.path.join(outdir, "trimmed_reads", os.path.basename(os.path.splitext(R2)[0]).replace(".fastq.gz", "").replace(".fq.gz", "") + ".trim.fastq.gz")
        spades_out=os.path.join(assembly_dir,self.sampleID+"_spades")

         # create metaspades command
        cmd= ["metaspades.py", "-1" , input_1, "-2", input_2,  "--only-assembler", "-t", "8", "-k", "21,31,41,51,61,71,81,91,101", "-o", spades_out]

        # log the spades command
        self.logger.info(f"Running command: {' '.join(cmd)}")
        print(f"Running metaspades on {self.sampleID}")
        if not os.path.isfile(os.path.join(spades_out, "contigs.fasta")):
            try:
                # run metaspades and capture output/error messages to logs
                subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
            except subprocess.CalledProcessError as error:
                # log errors and raise exception if subprocess fails
                self.logger.error(f"Error while running command: {' '.join(cmd)}")
                self.logger.error(f"Error message: {error.stderr}")
                raise
            print(f"Finished assembling sample {self.sampleID}")
            self.logger.info(f"Finished assembling sample {self.sampleID}")
        else:
            self.logger.info(f"Assembly for {self.sampleID} already exists. Skipping.")

    
    def blast(self, outdir, db, blastn):
        if blastn:
            method="BLASTn and BLAST"
        else:
            method="BLASTx and Diamond"
        
        print("\nSTEP 3: BLAST")
        # logger = logging.getLogger(self.sampleID)
        blast_dir=os.path.join(outdir,"blast_results")
        if not os.path.exists(blast_dir):
            os.mkdir(blast_dir)

        # log the start of the blast process
        print(f"Performing BLAST of sample {self.sampleID} with {method}, and saving files to {blast_dir}")
        self.logger.info(f"Database for BLAST located at {db}")
        self.logger.info(f"Performing BLAST of sample {self.sampleID} with {method}, and saving files to {blast_dir}")

        blast_query=os.path.join(outdir,"assemblies",self.sampleID+"_spades","contigs.fasta")
        blastn_out=os.path.join(blast_dir,self.sampleID+".blastn.results")
        blastx_out=os.path.join(blast_dir,self.sampleID+".blastx.results")

        # Create the blast commands
        blastn_cmd = ['blastn',
       '-query', blast_query,
       '-db', db,
       '-outfmt', '6 qseqid qlen sacc slen pident length evalue qstart qend sstart send sskingdoms salltitles',
       '-max_target_seqs', '1',
       '-max_hsps', '1',
       '-num_threads', '23',
       '-out', blastn_out]
        
        blastx_cmd = ["diamond",
        "blastx",
        "-d", db,
        "-q", blast_query,
        "--sensitive",
        "--outfmt", '6', 'qseqid', 'qlen', 'length', 'evalue', 'qstart', 'qend', 'sstart', 'send', 'stitle', 'salltitles',
        "--max-target-seqs", "1",
        "--max-hsps", "1",
        "--threads", "23",
        "-o", blastx_out]
        
        # Run the blastn command if blastn selected, otherwise run blastx
        if blastn:
            cmd=blastn_cmd
        else:
            cmd=blastx_cmd

        # log the blast command
        self.logger.info(f"Running command: {' '.join(cmd)}")
        print(f"Running command: {' '.join(cmd)}")
        print(f"Running {method} on {self.sampleID}")
        if blastn:
            if not os.path.isfile(blastn_out):
                try:
                    # run blastn and capture output/error messages to logs
                    subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
                except subprocess.CalledProcessError as error:
                    # log errors and raise exception if subprocess fails
                    self.logger.error(f"Error while running command: {' '.join(cmd)}")
                    self.logger.error(f"Error message: {error.stderr}")
                    raise
                print(f"Finished running BLASTn on sample {self.sampleID}")
                self.logger.info(f"Finished running BLASTn on sample {self.sampleID}")
            else:
                self.logger.info(f"BLASTn results for {self.sampleID} already exists. Skipping.")
        # Run BLASTx with diamond if BLASTn wasn't selected
        else:
            if not os.path.isfile(blastx_out):
                try:
                    # run blastx and capture output/error messages to logs
                    subprocess.run(cmd, check=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding='utf-8')
                except subprocess.CalledProcessError as error:
                    # log errors and raise exception if subprocess fails
                    self.logger.error(f"Error while running command: {' '.join(cmd)}")
                    self.logger.error(f"Error message: {error.stderr}")
                    raise
                print(f"Finished running BLASTx on sample {self.sampleID}")
                self.logger.info(f"Finished running BLASTx on sample {self.sampleID}")
            else:
                self.logger.info(f"BLASTx results for {self.sampleID} already exists. Skipping.")
            
    def run_pipeline(self):
        print(f"\n=== Starting Sample {self.sampleID} ===\n")
        self.trimming(self.R1, self.R2, self.adapters, self.threads, self.output)
        self.assembly(self.R1, self.R2, self.output)
        self.blast(self.output, self.db, self.blastn)
