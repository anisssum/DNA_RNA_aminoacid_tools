SHORT_CODE = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'U', 'O',
                'a', 'r', 'n', 'd', 'c', 'e', 'q', 'g', 'h', 'i', 'l', 'k', 'm', 'f', 'p', 's', 't', 'w', 'y', 'v', 'u', 'o']
LONG_CODE = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu', 'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu',
            'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', 'U': 'Sec', 'O': 'Pyl',
            'a': 'Ala', 'r': 'Arg', 'n': 'Asn', 'd': 'Asp', 'c': 'Cys', 'e': 'Glu', 'q': 'Gln', 'g': 'Gly', 'h': 'His', 'i': 'Ile', 'l': 'Leu',
            'k': 'Lys', 'm': 'Met', 'f': 'Phe', 'p': 'Pro', 's': 'Ser', 't': 'Thr', 'w': 'Trp', 'y': 'Tyr', 'v': 'Val', 'u': 'Sec', 'o': 'Pyl'}
MASS = {'A': 71.08, 'R': 156.2, 'N': 114.1, 'D': 115.1, 'C': 103.1, 'E': 129.1, 'Q': 128.1, 'G': 57.05, 'H': 137.1, 'I': 113.2, 'L': 113.2,
        'K': 128.2, 'M': 131.2, 'F': 147.2, 'P': 97.12, 'S': 87.08, 'T': 101.1, 'W': 186.2, 'Y': 163.2, 'V': 99.13, 'U': 168.05, 'O': 255.3,
        'a': 71.08, 'r': 156.2, 'n': 114.1, 'd': 115.1, 'c': 103.1, 'e': 129.1, 'q': 128.1, 'g': 57.05, 'h': 137.1, 'i': 113.2, 'l': 113.2,
        'k': 128.2, 'm': 131.2, 'f': 147.2, 'p': 97.12, 's': 87.08, 't': 101.1, 'w': 186.2, 'y': 163.2, 'v': 99.13, 'u': 168.05, 'o': 255.3}


def molecular_weight(seq: str) -> float:
    """
    Function calculates molecular weight of the amino acid chain
        Parameters:
            seq (str): each letter refers to one-letter coded proteinogenic amino acids
    Returns:
        (float) Molecular weight of tge given amino acid chain in Da
    """
    m = 0
    for acid in seq:
        m += MASS[acid]
    return m


def three_letter_code(seq: str) -> str:
    """
    Function converts single letter translations to three letter translations
        Parameters:
            seq (str): each letter refers to one-letter coded proteinogenic amino acids
        Returns:
            (str) translated in three-letter code
    """
    recording = seq.maketrans(LONG_CODE)
    return seq.translate(recording)


def show_length(seq: str) -> int:
    """
    Function counts the number of amino acids in the given sequence
        Parameters:
            seq (str): amino acid sequence
        Returns:
            (int): integer number of amino acid residues
    """
    return len(seq)


def folding(seq: str) -> str:
    """
    Counts the number of amino acids characteristic separately for alpha helixes and beta sheets,
    and gives out what will be the structure of the protein more.
    This function has been tested on proteins such as 2M3X, 6DT4 (PDB ID) and MHC, CRP.
    The obtained results corresponded to reality.
        Parameters:
            seq (str): amino acid sequence
        Returns:
            (str): overcoming structure ('alfa_helix', 'beta_sheet', 'equally')
    """
    alfa_helix = ['A', 'E', 'L', 'M', 'G', 'Y', 'S', 'a', 'e', 'l', 'm', 'g', 'y', 's']
    beta_sheet = ['Y', 'F', 'W', 'T', 'V', 'I', 'y', 'f', 'w', 't', 'v', 'i']
    alfa_helix_counts = 0
    beta_sheet_counts = 0
    for amino_acid in seq:
        if amino_acid in alfa_helix:
            alfa_helix_counts += 1
        elif amino_acid in beta_sheet:
            beta_sheet_counts += 1
    if alfa_helix_counts > beta_sheet_counts:
        return 'alfa_helix'
    elif alfa_helix_counts < beta_sheet_counts:
        return 'beta_sheet'
    elif alfa_helix_counts == beta_sheet_counts:
        return 'equally'


def seq_charge(seq: str) -> str:
    """
    Function evaluates the overall charge of the aminoacid chain in neutral aqueous solution (pH = 7)
        Parameters:
            seq (str): amino acid sequence of proteinogenic amino acids
        Returns:
            (str): "positive", "negative" or "neutral"
    Function realized by Anna Chesnokova
    """
    aminoacid_charge = {'R': 1, 'D': -1, 'E': -1, 'K': 1, 'O': 1, 'r': 1, 'd': -1, 'e': -1, 'k': 1, 'o': 1}
    charge = 0
    for aminoacid in seq:
        if aminoacid in aminoacid_charge:
            charge += aminoacid_charge[aminoacid]
    if charge > 0:
        return 'positive'
    elif charge < 0:
        return 'negative'
    else:
        return 'neutral'


def aminoacid_seqs_only(seqs: list) -> list:
    """
    Leaves only the amino acid sequences from the fed into the function.
        Parameters:
            seqs (list): amino acid sequence list
        Returns:
            aminoacid_seqs (list): amino acid sequence list without non amino acid sequence
    """
    aminoacid_seqs = []
    for seq in seqs:
        unique_chars = set(seq)
        amino_acid = set(SHORT_CODE)
        if unique_chars <= amino_acid:
            aminoacid_seqs.append(seq)
    return aminoacid_seqs
