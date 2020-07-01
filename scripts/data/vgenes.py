import numpy as np
import pandas as pd

# which are the correct columns?
ref_gene = 'IGHV1-46*01'
ref_cdr1 = 'GYTFTSYY'
ref_cdr2 = 'INPSGGST'
gap_character = '.'
gap_distance = 4
X_distance = 4 # penalty for aligning with an 'X' or with a '*' dunno what a '*' means maybe stop codon?

fasta_file = 'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP'

desired_species = 'Homo sapiens'
desired_prefix = 'IGHV'

# dictionary storing the aligned fasta sequences (amino acids)
alfas = {}

for line in open(fasta_file,'r'):
    if line[0] == '>':
        l = line.split('|')
        gene = l[1]
        species = l[2]
        if species!=desired_species or not gene.startswith(desired_prefix):
            gene = None
            continue
        alfas[gene] = ''

    elif gene:
        alfas[gene] += line.split()[0]

print(len(alfas), list(alfas.keys())[0])


alseq = alfas[ref_gene]
seq= alseq.replace(gap_character,'')
pos1 = seq.index(ref_cdr1)
pos2 = seq.index(ref_cdr2)


for i in range(len(alseq)):
    partseq = alseq[:i].replace(gap_character,'')
    if len(partseq) == pos1:
        cdr1_begin = i
    elif len(partseq) == pos1+len(ref_cdr1):
        cdr1_end = i
    elif len(partseq) == pos2:
        cdr2_begin = i
    elif len(partseq) == pos2+len(ref_cdr2):
        cdr2_end = i

print(cdr1_begin, cdr1_end, cdr2_begin, cdr2_end)

assert alfas[ref_gene][cdr1_begin:cdr1_end].replace(gap_character,'') == ref_cdr1
assert alfas[ref_gene][cdr2_begin:cdr2_end].replace(gap_character,'') == ref_cdr2

# read the tcrdist aa matrix
infofile = 'tcrdist_info_both_chains.txt'
dists = []
aas = []
for line in open(infofile,'r'):
    l = line.split()
    if l[0] == 'AAdist':
        aas.append(l[1])
        dists.append( [float(x) for x in l[2:] ] )
        assert len(dists[-1]) == 20
# convert to a dictionary indexed by aa pairs
aa_dist = {}

for i,a in enumerate(aas):
    for j,b in enumerate(aas):
        aa_dist[(a,b)] = dists[i][j]
        assert dists[i][j] == dists[j][i]
    aa_dist[('X',a)] = X_distance
    aa_dist[('*',a)] = X_distance
    aa_dist[(a,'X')] = X_distance
    aa_dist[(a,'*')] = X_distance
aa_dist[('X','*')] = X_distance
aa_dist[('*','X')] = X_distance

genes = sorted( alfas.keys())

cdr1s = [ alfas[x][cdr1_begin:cdr1_end] for x in genes]
cdr2s = [ alfas[x][cdr2_begin:cdr2_end] for x in genes]

for gene, cdr1, cdr2 in zip( genes, cdr1s, cdr2s):
    if 'X' in cdr1 or 'X' in cdr2:
        print('WHOAH:', gene, cdr1, cdr2)

def align_dist(a,b):
    assert len(a) == len(b)
    dist=0
    for aa,bb in zip(a,b):
        if aa==bb:
            continue
        elif aa==gap_character or bb==gap_character:
            dist += gap_distance
        else:
            dist += aa_dist[(aa,bb)]
    return dist

gene_dists = []
print('computing gene dists')
for ii in range(len(genes)):
    ii_dists = []
    for jj in range(len(genes)):
        dist = align_dist(cdr1s[ii], cdr1s[jj]) + align_dist(cdr2s[ii], cdr2s[jj])
        ii_dists.append(dist)
    gene_dists.append(ii_dists)
gene_dists = np.array(gene_dists)
print(gene_dists.shape)

df = pd.DataFrame( { g:gene_dists[ii,:] for ii,g in enumerate(genes) } )
df['gene'] = genes
df.set_index('gene', inplace=True)

print(df.head())

df['cdr1'] = cdr1s
df['cdr2'] = cdr2s

df = df[ ['cdr1','cdr2']+genes ]

df.to_csv('all_IGHV_dists.tsv', sep='\t')


