from Data.DataUtils.Utils import read_fasta
import os

'''
adding M as first amino acid where necessary. 
'''

# read in sampled sequences after identity filtering and clustering
sampled_seqs_path = '../Data/samples_filterstep2.fasta'
seqs = read_fasta(sampled_seqs_path)
altered = {}
final_seqs = {}

# add M where necessary
for name, seq in seqs.items():
    if seq[0] != 'M':
        altered[name] = 'M' + seq
        final_seqs[name] = 'M' + seq
    else:
        final_seqs[name] = seq

# save altered sequences to fasta file "final_sequences.fasta"
content = '\n'.join(['>' + head + '\n' + seq for head, seq in altered.items()])
open('../Data/final_samples.fasta', 'w').write(content)

# save altered seqs to fasta files, one dimer and one MinD-interaction
af_dir_path = '../Data/AlphaFold'
if not os.path.isdir(af_dir_path):
    os.mkdir(af_dir_path)
minD = read_fasta('./Data/EColi_MinD.fasta')
# MinD-interaction
for head, seq in altered.items():
    content = '>' + head + '\n' + seq + '\n>' + list(minD.keys())[0] + '\n' + list(minD.values())[0]
    open(os.path.join(af_dir_path, head + '.fasta'), 'w').write(content)
# dimer
for head, seq in altered.items():
    content = '>' + head + '\n' + seq
    content += '\n' + content
    open(os.path.join(af_dir_path, head + '_dimer.fasta'), 'w').write(content)

