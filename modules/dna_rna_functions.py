DNA_COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 't', 't': 'a', 'c': 'g', 'g': 'c'}
RNA_COMPLEMENT = {'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C', 'a': 'u', 'u': 'a', 'c': 'g', 'g': 'c'}


def dna_or_rna_only(seqs):
    """
        The function leaves only DNA and RNA sequences
            Parameters:
                seqs - sequences fed to the input of the main function
            Return:
                list of DNA and RNA sequences only
    """
    dna_seqs = []
    for seq in seqs:
        unique_chars = set(seq)
        nuc = set('ATGCatgcUu')
        if unique_chars <= nuc:
            dna_seqs.append(seq)
    return dna_seqs


def transcribe(seq):
    """
        The function returns the transcribed sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) transcript sequence
    """
    transcript = seq.replace('T', 'U').replace('t', 'u')
    return transcript


def reverse(seq):
    """
        The function produces a reverse sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) reverse sequence
    """
    reverse_seq = seq[::-1]
    return reverse_seq


def complement(seq):
    """
        The function produces a complementary sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) complementary sequence
    """
    complement_seq = []
    if 'U' in set(seq):
        using_dict = RNA_COMPLEMENT
    else:
        using_dict = DNA_COMPLEMENT
    for nucl in seq:
        complement_seq.append(using_dict[nucl])
    return ''.join(complement_seq)


def reverse_complement(seq):
    """
        The function produces a reverse and complementary sequence
            Parameters:
                seq (str) - sequence
            Return:
                (str) reverse and complementary sequence
    """
    reverse_compl_seq = complement(reverse(seq))
    return reverse_compl_seq
