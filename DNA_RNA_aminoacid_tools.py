import modules.fastq_filtering as ff
import modules.aminoacid_functions as af
import modules.dna_rna_functions as drf

def filter_seqs(input_path: str, gc_bounds: tuple, length_bounds: tuple, quality_threshold = 0, output_filename = ''):
    """
    Performs functions for working with fastq.
        Parameters:
            - input_path - takes as input the path to the FASTQ-file.
                        A DNA or RNA sequences can consist of either uppercase or lowercase letters.
            - gc_bounds - GC composition interval (in percent) for filtering, by default it is (0, 100).
                        If you pass one number as an argument, it is considered to be the upper limit.
                        Examples: gc_bounds = (20, 80) - save only reads with GC content from 20 to 80%,
                        gc_bounds = 44.4 - save reads with GC content less than 44.4%.
            - length_bounds - length interval for filtering, everything is similar to gc_bounds,
                        but by default it is equal to (0, 2**32).
            - quality_threshold - threshold value of average read quality for filtering, default is 0 (phred33 scale).
                        Reads with average quality across all nucleotides below the threshold are discarded.
            - output_filename - name of the file with the filtering result, if not specified, the name of the input file is assigned by default.
                        Important: specify the file extension.
        Example input:
            filter_seqs(seqs = 'example_fastq.fastq', gc_bounds = (20, 80), length_bounds = (0, 89), quality_threshold = 34), output_filename = 'my_out_file'
        Return:
            The function returns a file consisting of only those sequences that satisfy all conditions.
            It puts this file into the "fastq_filtrator_resuls" folder and creates it, if it does not exist.
            All described intervals include both upper and lower boundaries.
            Depending on the function being performed, the following returns will occur:
                If the sequences in the dictionary are RNA, then there will be no filtering by gc composition.
                If you supply a tuple of more than 2 values for the gc_bounds and length_bounds arguments,
                you will receive the errors "Incorrect gc_bounds input" and "Incorrect length_bounds input" respectively.
    """
    filtered_dict = ff.read_fastq(input_path)
    if type(gc_bounds) is tuple and len(gc_bounds) == 2:
        min_gc_bound = gc_bounds[0]
        max_gc_bound = gc_bounds[1]
    elif (type(gc_bounds) is int) or (type(gc_bounds) is float):
        min_gc_bound = 0
        max_gc_bound = gc_bounds
    else:
        return print("Incorrect gc_bounds input")
    if type(length_bounds) is tuple and len(length_bounds) == 2:
        min_len_bound = length_bounds[0]
        max_len_bound = length_bounds[1]
    elif type(length_bounds) is tuple and len(length_bounds) == 0:
        min_len_bound = 0
        max_len_bound = 2**32
    elif (type(length_bounds) is int) or (type(length_bounds) is float):
        min_len_bound = 0
        max_len_bound = length_bounds
    else:
        return print("Incorrect length_bounds input")
    func_dict = {ff.filter_gc: (min_gc_bound, max_gc_bound),
                 ff.filter_length: (min_len_bound, max_len_bound),
                 ff.filter_quality_threshold: (quality_threshold,)}
    for func in func_dict:
        key_list = []
        for filtered_dict_key in filtered_dict:
            key_list.append(filtered_dict_key)
        for key in key_list:
            if len(func_dict[func]) == 2:
                func(filtered_dict, key, func_dict[func][0], func_dict[func][1])
            else:
                func(filtered_dict, key, func_dict[func][0])
    if output_filename == '':
        output_filename = input_path.split('/')[-1]
    return ff.save_fastq(filtered_dict, output_filename)


def amino_acid_tools(*args: str):
    """
    Performs functions for working with protein sequences.
        Parameters:
            The function must accept an unlimited number of protein sequences (str) as input,
            the last  variable must be the function (str) you want to execute.
            The amino acid sequence can consist of both uppercase and lowercase letters.
        Input example:
            amino_acid_tools('PLPKVEL','VDviRIkLQ','PPDFGKT','folding')
        Functions:
            molecular_weight: calculates molecular weight of the amino acid chain
            three_letter_code: converts single letter translations to three letter translations
            show_length: counts the number of amino acids in the given sequence
            folding: counts the number of amino acids characteristic separately for alpha helixes and beta sheets,
                    and gives out what will be the structure of the protein more
            seq_charge: evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7)
        Returns:
            If one sequence is supplied, a string with the result is returned.
            If several are submitted, a list of strings is returned.
            Depending on the function performed, the following returns will occur:
                molecular_weight (int) or (list): amino acid sequence molecular weight number or list of numbers
                three_letter_code (str) or (list): translated sequence from one-letter in three-letter code
                show_length (int) or (list): integer number of amino acid residues
                folding (str) or (list): 'alpha_helix', if there are more alpha helices
                                        'beta_sheet', if there are more beta sheets
                                        'equally', if the probability of alpha spirals and beta sheets are the same
                			seq_charge(str) or (list): "positive", "negative" or "neutral"
    """
    *seqs, function = args
    d_of_functions = {'molecular_weight': af.molecular_weight,
                      'three_letter_code': af.three_letter_code,
                      'show_length': af.show_length,
                      'folding': af.folding,
                      'seq_charge': af.seq_charge}
    answer = []
    aminoacid_seqs = af.aminoacid_seqs_only(seqs)
    for sequence in aminoacid_seqs:
        answer.append(d_of_functions[function](sequence))
    if len(answer) == 1:
        return answer[0]
    return answer


def run_dna_rna_tools(*args: str):
    """
        Performs functions for processing one or several DNA and/or RNA sequences.
           Parameters:
               The function must accept an unlimited number of sequences (str) as input.
               the last variable should be the function (str) you want to execute.
               The sequence can consist of both uppercase and lowercase letters.
           Example input:
               def run_dna_rna_tools("ATgAAaC", "cUgAuaC", "reverse")
           Functions:
               transcribe — print the transcribed sequence
               reverse — print the reversed sequence
               complement — print the complementary sequence
               reverse_complement — print the reverse complementary sequence
           Return:
               If a single sequence is specified, a string containing the result is returned.
               If multiple strings are sent, a list of strings is returned.
               Depending on the function being performed, the following returns will occur:
                   - if the sequences you pass are not DNA or RNA, then the result of the function will be a list without them
                   - if the sequence is RNA, then the 'transcribe' procedure will produce the unchanged sequence
    """
    *seqs, function = args
    d_of_functions = {'transcribe': drf.transcribe, 
                      'reverse': drf.reverse,
                      'reverse_complement': drf.reverse_complement,
                      'complement': drf.complement,
                     }
    answer = []
    dna_rna_seqs = drf.dna_or_rna_only(seqs)
    for sequence in dna_rna_seqs:
        answer.append(d_of_functions[function](sequence))
    if len(answer) == 1:
        return answer[0]
    return answer
