# Amino acid, DNA and RNA sequence processing and fastq sequences filtering

This tool is designed to work with fastq, fasta ang gbk files, amino acid sequences, DNA and RNA sequences.

## Usage

Depending on your data, you can use three functions from `DNA_RNA_aminoacid_tools.py` main function: `filter_seqs`, `amino_acid_tools`, `run_dna_rna_tools`.

And two functions from `bio_files_processor.py` main function: `convert_multiline_fasta_to_oneline` and `select_genes_from_gbk_to_fasta`.

If you are working with fastq file, you call the `filter_seqs` function, that takes as input the fastq file containing an unlimited number of fastq sequences, filtering bounds on gz composition, length, and average read quality. The function then filters all sequences in the file. Sequences can contain both uppercase and lowercase letters. The result of the function is the output file with filtered sequences.

If you are working with amino acid sequences, you call the `amino_acid_tools` function, which takes as input an arbitrary number of arguments with amino acid sequences (str) as well as the name of the procedure to be performed (it is always the last argument, str). The command then performs the specified action on all the given sequences. If a single sequence is submitted, a string with the result is returned. If multiple sequences are supplied, a list of strings is returned.

If you are working with dna or rna sequences, you call `run_dna_rna_tools` takes as input an arbitrary number of arguments with DNA or RNA sequences (str) and the name of the procedure to be performed (it is always the last argument, str). After that the command performs the specified action on all the given sequences. If one sequence is submitted, a string with the result is returned. If several sequences are submitted, a list of strings is returned.

For `amino_acid_tools` and `run_dna_rna_tools` input sequences can contain both uppercase and lowercase letters, but the last argument with the function name must match the possible functions.

If you are working with fasta file, you call the `convert_multiline_fasta_to_oneline` function, that takes as input the fasta file containing an unlimited number of fasta sequences. In which the sequence (DNA/RNA/protein/ ... ) can be split into several lines, and then saves it into a new fasta file in which each sequence fits into a single line.

If you are working with gdk file, you call the `select_genes_from_gbk_to_fasta` function, that takes as input the gbk file. The function selects some number of genes before and after each of the genes of interest and stores their protein sequence (translation) into a fasta file.

### Remark

#### For `filter_seqs` function

- if the sequences in the file are RNA, then there will be no filtering by gc composition
- if you supply a tuple of more than 2 values for the gc_bounds and length_bounds arguments, you will receive the errors "Incorrect gc_bounds input" and "Incorrect length_bounds input" respectively

#### For `amino_acid_tools` function
- if the sequences passed by you contain inappropriate characters (not from the single-letter aminoxylot encoding), the result of the function will be a list without them
- the fewer amino acids a sequence contains, the less reliable the 'folding' function is

#### For `run_dna_rna_tools` function
- if the sequences you passed are not DNA or RNA, the result of the function will be a list without them
- if the sequence is RNA, the 'transcribe' procedure will produce an unmodified sequence

#### For `select_genes_from_gbk_to_fasta` function
- If the gene of interest is extreme, there will be no environment on one side

## Parameters 

#### For `filter_seqs` function

- **input_path**: takes as input the path to the FASTQ-file. A DNA or RNA sequences can consist of either uppercase or lowercase letters. Examples: input_path = 'example_fastq.fastq'.
- **gc_bounds**: GC composition interval (in percent) for filtering, by default it is (0, 100). If you pass one number as an argument, it is considered to be the upper limit. Examples: gc_bounds = (20, 80): save only reads with GC content from 20 to 80%, gc_bounds = 44.4: save reads with GC content less than 44.4%.
- **length_bounds**: length interval for filtering, everything is similar to gc_bounds, but by default it is equal to (0, 2**32).
- **quality_threshold**: threshold value of average read quality for filtering, default is 0 (phred33 scale). Reads with average quality across all nucleotides below the threshold are discarded. 
- **output_filename**: name of the file with the filtering result, if not specified, the name of the input file is assigned by default.

#### For `amino_acid_tools` function

- **molecular_weight**: calculate the molecular weight of the amino acid chain in Da, according to the average amino acid residues molecular masses rounded to 1 or 2 decimal places.
- **three_letter_code**: converts standard single letter translations to three letter translations.
- **show_length**: count the overall number of amino acids in the given.
- **sequence folding**: count the number of amino acids characteristic separately for alpha helixes and beta sheets,and give out what will be the structure of the protein more. This function has been tested on proteins such as 2M3X, 6DT4 (PDB ID) and MHC, CRP. The obtained results corresponded to reality.
- **seq_charge**: evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7), according to the pKa of amino acid side chains, lysine, pyrrolizine and arginine contribute +1, while asparagine and glutamic amino acids contribute -1. The total charge of a protein is evaluated as positive, negative, or neutral as the sum of these contributions.

#### For `run_dna_rna_tools` function

- **transcribe**: print transcribed sequence
- **reverse**: print reversed sequence
- **complement**: print complementary sequence
- **reverse_complement**: print reverse complementary sequence

#### For `convert_multiline_fasta_to_oneline` function

- **input_fasta**: the path to the FASTA file is taken as input
- **output_fasta**: name of output file, if not specified then adds out to the name of the source file at the beginning

#### For `select_genes_from_gbk_to_fasta' function

- **input_gbk** - the path to GBK file is taken as input data
- **genes** - genes of interest, near which neighbors are searched for
- **n_before**, **n_after** - number of genes before and after (>0), default values - 1
- **output_fasta** - the name of the output file, if not specified, out is appended to the name of the source file at the beginning


## Examples

### For `filter_seqs` function

Below is an example of processing a fastq sequences.
```shell
EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079804:1:SRR292678:1:1101:175500:175500': ('GACGCCGTGGCTGCACTATTTGAGGCACCTGTCCTCGAAGGGAAGTTCATCTCGACGCGTGTCACTATGACATGAATG', 'GGGGGFFCFEEEFFDGFBGGGA5DG@5DDCBDDE=GFADDFF5BE49<<<BDD?CE<A<8:59;@C.C9CECBAC=DE'),
    '@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT', 'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE'),
    '@SRX079804:1:SRR292678:1:1101:190845:190845': ('CCTCAGCGTGGATTGCCGCTCATGCAGGAGCAGATAATCCCTTCGCCATCCCATTAAGCGCCGTTGTCGGTATTCC', 'FF@FFCFEECEBEC@@BBBBDFBBFFDFFEFFEB8FFFFFFFFEFCEB/>BBA@AFFFEEEEECE;ACD@DBBEEE'),
    '@SRX079804:1:SRR292678:1:1101:198993:198993': ('AGTTATTTATGCATCATTCTCATGTATGAGCCAACAAGATAGTACAAGTTTTATTGCTATGAGTTCAGTACAACA', '<<<=;@B??@<>@><48876EADEG6B<A@*;398@.=BB<7:>.BB@.?+98204<:<>@?A=@EFEFFFEEFB'),
    '@SRX079804:1:SRR292678:1:1101:204480:204480': ('AGTGAGACACCCCTGAACATTCCTAGTAAGACATCTTTGAATATTACTAGTTAGCCACACTTTAAAATGACCCG', '<98;<@@@:@CD@BCCDD=DBBCEBBAAA@9???@BCDBCGF=GEGDFGDBEEEEEFFFF=EDEE=DCD@@BBC')
    }
filter_seqs(seqs = EXAMPLE_FASTQ, gc_bounds = (20, 80), length_bounds = (0, 89), quality_threshold = 34)
```

Input: EXAMPLE_FASTQ, gc_bounds = (20, 80), length_bounds = (0, 89), quality_threshold = 34

Output: {'@SRX079804:1:SRR292678:1:1101:190136:190136': ('GAACCTTCTTTAATTTATCTAGAGCCCAAATTTTAGTCAATCTATCAACTAAAATACCTACTGCTACTACAAGTATT',
'DACD@BEECEDE.BEDDDDD,>:@>EEBEEHEFEHHFFHH?FGBGFBBD77B;;C?FFFFGGFED.BBABBG@DBBE')}

### For `amino_acid_tools` function

Below is an example of processing an amino acid sequence.

#### Using the function for molecular weight calculation

```shell  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'molecular_weight')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'molecular_weight'

Output: [1228.66, 1447.8400000000001, 1224.6399999999999]

#### Using the function to convert one-letter translations to three-letter translations

```shell  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'three_letter_code')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'three_letter_code'

Output: ['GluGlyValIleMetSerGluLeuLysLeuLys', 'ProLeuProLysValGluLeuProProAspPheValAsp', 'AspValIleGlyIleSerIleLeuGlyLysGluVal']

#### Using the function to counts the number of amino acids in the given sequence

```shell  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'show_length')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'show_length'

Output: [11, 13, 12]

#### Using the function to determine the predominant secondary structure

```shell  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'folding')  
```
Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'folding'

Output: ['alfa_helix', 'equally', 'equally']

#### Using the function to estimate relative charge

```shell  
amino_acid_tools('EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'seq_charge')  
```

Input: 'EGVIMSELKLK', 'PLPKvelPPDFVD', 'DVIGISILGKEV', 'seq_charge'

Output: ['neutral', 'negative', 'negative']

#### For `run_dna_rna_tools` function

Below is an example of processing DNA and RNA sequences.

#### Using the function to reverse sequence

```shell 
run_dna_rna_tools('ATGCccT', 'UGACU', 'reverse')
```
Input: 'ATGCccT', 'UGACT', 'reverse')

Output: ['TccCGTA', 'UCAGU']

#### Using the function to transcribe sequence

```shell 
run_dna_rna_tools('ATGCccT', 'UGACU', 'reverse')
```
Input: 'ATGCccT', 'UGACT', 'reverse')

Output: ['AUGCccU', 'UGACU']

#### Using the function to complement sequence

```shell 
run_dna_rna_tools('ATGCccT', 'UGACU', 'reverse')
```
Input: 'ATGCccT', 'UGACT', 'complement')

Output: ['TACGggA', 'ACUGA']

#### Using the function to reverse and complement sequence

```shell 
run_dna_rna_tools('ATGCccT', 'UGACU', 'reverse')
```
Input: 'ATGCccT', 'UGACT', 'reverse')

Output: ['AggGCAT', 'AGUCA']

#### For `convert_multiline_fasta_to_oneline` function

```shell 
convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta', output_fasta = 'qwe')
```
Output: 'out_qwe.fasta' file

#### For `select_genes_from_gbk_to_fasta` function

```shell 
select_genes_from_gbk_to_fasta('example_gbk.gbk', 'phrB', 'galK', n_before = 1, n_after = 1, output_fasta = '')
```
Output: 'out_example_gbk.fasta' file

## Requirements

You need to import modules for 'DNA_RNA_aminoacid_tools.py':

- import modules.fastq_filtering as ff
- import modules.aminoacid_functions as af
- import modules.dna_rna_functions as drf

For `bio_files_processor.py`:

- import os

## Contacts  
- [Cesnokova Anna] annachesnokova0303@gmail.com
