# Homework for the spring semester of the Bioinformatics Institute

This tool contains the following topics: OOP, iterators, decorators, API, parallel programming.

## Installation

1. Clone this repository:

```shell
git clone git@github.com:anisssum/HW18
```

2. Navigate to the cloned repository:

```shell
cd HW18
```

3. Install the required dependencies:

```shell
pip install -r requirements.txt
```

## Tools

### DNA_RNA_aminoacid_tools

This script provides various functionalities for processing biological sequences, including DNA, RNA, and protein sequences. It also offers tools for analyzing and manipulating sequences, as well as performing tasks such as filtering FASTQ files and running GENSCAN analysis.

### bio_files_processor

This script provides functionalities to work with FASTA files containing DNA, RNA, and protein sequences. It includes the following functions:

- convert_multiline_fasta_to_oneline: Converts a FASTA file where each sequence spans multiple lines into a file where each sequence is on a single line.

- select_genes_from_gbk_to_fasta: Selects a specified number of genes before and after each gene of interest from a GenBank (GBK) file and stores their protein sequences into a FASTA file.

- FastaRecord class: Defines a data class representing a FASTA record, containing an ID, description, and sequence.

- OpenFasta class: Implements a context manager to open a FASTA file, iterate through its records, and read records.

### custom_random_forest

This script implements a custom RandomForestClassifier, which is an ensemble learning method for classification that works by building multiple decision trees during training and outputs a class that is the mode of the classes (classification) or the average prediction (regression) of the individual trees. Its feature is that you can parallele the processes by selecting the desired number of kernels.

### Showcases

Contains examples of each tool's work.

### test_biotools

This script contains unit tests and assertions for functions related to DNA, RNA, and protein sequence manipulation, as well as file processing utilities.

## Author

- [Cesnokova Anna] annachesnokova0303@gmail.com
