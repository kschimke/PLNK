# PLNK
PLNK runs alongside an Oxford Nanopore MinION sequencer, processing individual fast5 files using guppy for basecalling, C3POa for R2C2 consensus calling, and mappy for alignment before analyzing the library content.

Expects R2C2 reads using samples previously prepped for Illumina sequencing. 

--------------------------------------------------------------------------------

## Dependencies

- [Python 3](https://www.python.org/downloads/)
- [NumPy](https://pypi.org/project/numpy/)
- [guppy](https://community.nanoporetech.com/downloads) (requires login)
- [C3Poa](https://github.com/christopher-vollmers/C3POa)
- [mappy](https://pypi.org/project/mappy/)

--------------------------------------------------------------------------------

## Usage

After resolving all of the dependencies, PLNK will take <3 minutes to install and runs with python.

## PLNK.py

```bash
python3 PLNK.py -i /path/to/fast5/directory/ 
                -o /path/to/output/directory/ 
                -s splint.fasta 
                -sm samplesheet.tsv
                -t target.bed
                -r reference.fasta
                -a adapter.fasta
                -c config 
                -C /path/to/C3Poa/directory/
                -n 1
                -g 0 
```

Arguments:
```
-i  directory collecting fast5 files

-o  output path

-s  sequence of DNA splint used in R2C2 protocol in fasta format

-sm sample sheet with entries for each sample in sequencing data

-t  regions of interest within the genome in bed format

-r  reference genome in fasta format

-a  sequence of cDNA adapter sequences in fasta format. Sequence names must be
    3Prime_adapter and 5Prime_adapter

-c  config file containing path to BLAT and racon binaries

-C  path to C3Poa directory

-n  number of threads to use (defaults to 1)

-g  cuda numbers for GPUs, comma separated list

-T  print time information for processing fast5 files

-V  verbose, prints tool output to console otherwise creates log files

-v  print the PLNK version and exit
```

Example output directory tree:
```
output_dir
├── Fast5_# 
│   ├── Consensus
│   │   ├── Demultiplexed
│   │   │   ├── Sample_Name.fasta
│   │   │   └── Undetermined.fasta
│   │   ├── Index#
│   │   │   ├── R2C2_full_length_consensus_reads.fasta
│   │   │   ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   │   │   ├── R2C2_full_length_consensus_reads_right_splint.fasta
│   │   │   ├── R2C2_Subreads.fastq
│   │   │   └── R2C2_Consensus.fasta
│   │   └── no_index_found
│   │       ├── R2C2_full_length_consensus_reads.fasta
│   │       ├── R2C2_full_length_consensus_reads_left_splint.fasta
│   │       └── R2C2_full_length_consensus_reads_right_splint.fasta
```
--------------------------------------------------------------------------------

## sequencing_simulator.py

Simulates a MinION sequencing run by copying fast5 files from the input directory to an 
output directory at a specified rate or by comparing the metadata of the input files.

```bash
python3 sequencing_simulator.py -i /path/to/fast5/directory/ 
                                -o /path/to/output/directory/ 
                                -r 1.5  
```
Arguments:
```
-i  directory containing fast5 files from a previous sequencing run

-o  directory fast5 files will be copied over to

-r  a float specifying the number of minutes between file transfers, 
    defaults to extrapolating a rate from the fast5s' metadata
```
