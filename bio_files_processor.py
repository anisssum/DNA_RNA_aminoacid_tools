import os


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta = ''):
    """
    Performs functions to work with fasta file with DNA, RNA and protein sequences.
        Parameters:
            - input_fasta - the path to the FASTA file is taken as input
            - output_fasta - name of output file, if not specified then adds out to the name of the source file at the beginning
        Example input:
            convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta', output_fasta = 'qwe')
        Return:
            The function returns saves the input file to a new fasta file in which each sequence fits on a single line.
    """
    d = {}
    if output_fasta == '':
        output_fasta = input_fasta.split('/')[-1]
    with open(input_fasta,'r') as fasta_file:
        fasta_list = fasta_file.readlines()
    with open('out_'+output_fasta+'.fasta','w') as output_fasta:
        i=0
        while i < len(fasta_list)-1:
            if fasta_list[i][0] == '>':
                output_fasta.write(fasta_list[i])
                i+=1
            else:
                y=i
                string = ''
                while y < len(fasta_list) and fasta_list[y][0] != '>':
                    string = string+fasta_list[y].replace('\n', '')
                    y+=1
                output_fasta.write(string+'\n')
                i=y


def select_genes_from_gbk_to_fasta(input_gbk: str, *genes: str, n_before = 1, n_after = 1, output_fasta = ''):
    """
    The function selects some number of genes before and after each of the genes of interest and stores their protein sequence (translation) into a fasta file.
        Parameters:
            - input_gbk - the path to GBK file is taken as input data
            - genes - genes of interest, near which neighbors are searched for
            - n_before, n_after - number of genes before and after (>0), default values - 1
            - output_fasta - the name of the output file, if not specified, out is appended to the name of the source file at the beginning
        Input example:
            select_genes_from_gbk_to_fasta('example_gbk.gbk', 'phrB', 'galK', n_before = 1, n_after = 1, output_fasta = '')
        Return:
            The function saves the sequences found to a new fasta file.
            If the gene of interest is extreme, there will be no environment on one side.
    """
    with open(input_gbk,'r') as gbk_file:
        gbk_list = gbk_file.readlines()
    all_genes = [s.split('"')[1] for s in gbk_list if "/gene=" in s]
    new_genes = []
    for gene in genes:
        ind_gene = all_genes.index(gene)
        i=1
        j=1
        while i <= n_before:
            if ind_gene-i >= 0:
                new_genes.append(all_genes[ind_gene-i])
            i+=1
        while j <= n_after:
            if ind_gene+j <= len(all_genes):
                new_genes.append(all_genes[ind_gene+j])
            j+=1
    indexes = []
    for new_gene in new_genes:
        indexes.append([gbk_list.index(s) for s in gbk_list if new_gene in s])
    seq = ''
    seqs = []
    for ind in indexes:
        start = ind[0]
        while '/translation=' not in gbk_list[start]:
            start+=1
        end = start+1
        seq = gbk_list[start].split('"')[1].replace('\n', '')
        while '"' not in gbk_list[end]:
            seq = seq+gbk_list[end].replace('\n', '').replace(' ', '')
            end+=1
        seq = seq+gbk_list[end].split('"')[0].replace(' ', '')
        seqs.append(seq)
    if output_fasta == '':
        output_fasta = input_gbk.split('/')[-1]
    with open('out_'+output_fasta+'.fasta','w') as output_fasta:
        for n in range(len(new_genes)):
            output_fasta.write('>'+new_genes[n]+'\n')
            output_fasta.write(seqs[n]+'\n')
