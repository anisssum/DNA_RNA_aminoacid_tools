# Amino acid, DNA and RNA sequence processing and fastq sequences filtering 
This tool is designed to work with fastq sequences, amino acid sequences, DNA and RNA sequences.

## Usage  
Depending on your data, you can use three functions: `filter_seqs`, `amino_acid_tools`, `run_dna_rna_tools`.

If you are working with fastq sequences, you call the `filter_seqs` function, that takes as input a dictionary containing an unlimited number of fastq sequences, filtering bounds on gz composition, length, and average read quality. The function then filters all sequences in the dictionary. Sequences can contain both uppercase and lowercase letters. The result of the function is a dictionary with filtered sequences.

If you are working with amino acid sequences, you call the `amino_acid_tools` function, which takes as input an arbitrary number of arguments with amino acid sequences (str) as well as the name of the procedure to be performed (it is always the last argument, str). The command then performs the specified action on all the given sequences. If a single sequence is submitted, a string with the result is returned. If multiple sequences are supplied, a list of strings is returned.

If you are working with dna or rna sequences, you call `run_dna_rna_tools` takes as input an arbitrary number of arguments with DNA or RNA sequences (str) and the name of the procedure to be performed (it is always the last argument, str). After that the command performs the specified action on all the given sequences. If one sequence is submitted, a string with the result is returned. If several sequences are submitted, a list of strings is returned.

For `amino_acid_tools` and `run_dna_rna_tools` input sequences can contain both uppercase and lowercase letters, but the last argument with the function name must match the possible functions.

### Remark

#### For `filter_seqs` function
- if the sequences in the dictionary are RNA, then there will be no filtering by gc composition.
- if you supply a tuple of more than 2 values for the gc_bounds and length_bounds arguments, you will receive the errors "Incorrect gc_bounds input" and "Incorrect length_bounds input" respectively.

#### For `amino_acid_tools` function
- if the sequences passed by you contain inappropriate characters (not from the single-letter aminoxylot encoding), the result of the function will be a list without them
- the fewer amino acids a sequence contains, the less reliable the 'folding' function is

#### For `run_dna_rna_tools` function
- if the sequences you passed are not DNA or RNA, the result of the function will be a list without them.
- if the sequence is RNA, the 'transcribe' procedure will produce an unmodified sequence

## Parameters 

#### For `filter_seqs` function
- **seqs**: a dictionary with an unlimited number of fastq sequences. The structure is as follows. Key - string, sequence name. The value is a tuple of two strings: sequence and quality. A DNA sequence can consist of either uppercase or lowercase letters.
- **gc_bounds**: GC composition interval (in percent) for filtering, by default it is (0, 100). If you pass one number as an argument, it is considered to be the upper limit. Examples: gc_bounds = (20, 80): save only reads with GC content from 20 to 80%, gc_bounds = 44.4: save reads with GC content less than 44.4%.
- **length_bounds**: length interval for filtering, everything is similar to gc_bounds, but by default it is equal to (0, 2**32).
- **quality_threshold**: threshold value of average read quality for filtering, default is 0 (phred33 scale). Reads with average quality across all nucleotides below the threshold are discarded.

#### For `amino_acid_tools` function
The following options for aminoacid sequence processing are available at the moment:

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

## Requirements
You need to import modules:

- import modules.fastq_filtering as ff
- import modules.aminoacid_functions as af
- import modules.dna_rna_functions as drf

## Contacts  
- [Cesnokova Anna] annachesnokova0303@gmail.com
