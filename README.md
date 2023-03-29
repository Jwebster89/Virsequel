# Virsequel

Virsequel is a Python script for identifying viral sequences in RNA-seq data. It uses a combination of adapter trimming and assembly to identify and extract viral sequences from raw RNA-seq reads.

## Requirements

To run Virsequel, you'll need the following software and data:

- Python 3
- BBMap
- MetaSPAdes
- A database of viral reference sequences (To be included)

## Installation

1. Clone the Virsequel repository from GitHub:

```
git clone https://github.com/Jwebster89/virsequel.git
```


2. Install the required software:

- Python 3: Follow the instructions for your operating system to install Python 3.
- BBMap: Download and install BBMap from the [BBMap website](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/).
- MetaSPAdes: Download and install MetaSPAdes from the [MetaSPAdes website](https://github.com/ablab/spades).

3. (To be included) Download a database of viral reference sequences to improve viral sequence identification. This step is optional, blastn can be performed with the "--remote" flag.


## Usage

To use Virsequel, run the `run_virsequel.py` script with the following command-line arguments:

```
usage: run_virsequel.py -1 READ_1 -2 READ_2 -a ADAPTER -o OUTPUT [-h] [-t THREADS]

        Identify viral sequences in rna seq data 

    Version: 0.1.0


Required Arguments:
  -1 READ_1, --read_1 READ_1
                        Forward RNA fastq.gz data
  -2 READ_2, --read_2 READ_2
                        Reverse RNA fastq.gz data
  -a ADAPTER, --adapter ADAPTER
                        adapter file for bbduk trimming
  -o OUTPUT, --output OUTPUT
                        output path

Optional Arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads for non-spades processes
```

Virsequel will generate trimmed reads and meta RNA assemblies in the specified output directory.

## License

Virsequel is released under the [GNU General Public License v3.0](https://opensource.org/licenses/GPL-3.0). See the `LICENSE` file for more information.

## Contributing

If you find a bug or have a feature request, please create an issue on the GitHub repository. If you'd like to contribute code, please fork the repository and create a pull request with your changes.

## Credits

Virsequel was created by [John Webster](https://github.com/Jwebster89).


