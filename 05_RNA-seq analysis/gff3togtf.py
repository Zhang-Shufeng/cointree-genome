import sys

inFile = open(sys.argv[1], 'r')

ID_gene = ''

for line in inFile:
    # skip comment lines that start with the '#' character
    if line[0] != '#':
        # split line into columns by tab
        data = line.strip().split('\t')

        ID_gene = ID_gene
        ID_mRNA = ''

        # if the feature is a gene
        if data[2] == "gene":
            # get the id
            ID_gene = data[-1].split('ID=')[-1].split(';')[0]
            data[-1] = 'gene_id "' + ID_gene + '"; '

        elif data[2] == "mRNA":
            # get two id
            ID_mRNA = data[-1].split('ID=')[-1].split(';')[0]
            ID_gene = data[-1].split('Parent=')[-1].split(';')[0]
            data[-1] = 'gene_id "' + ID_gene + '"; transcript_id "' + ID_mRNA + '"; '

        # if the feature is anything else
        else:
            # get the parent as the ID
            ID_mRNA = data[-1].split('Parent=')[-1].split(';')[0]
            data[-1] = 'gene_id "' + ID_gene + '"; transcript_id "' + ID_mRNA + '"; '

        # print out this new GTF line
        print('\t'.join(data))
