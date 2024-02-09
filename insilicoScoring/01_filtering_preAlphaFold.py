from Filtering.Filters import Filter
from Data.DataUtils.Utils import read_fasta
import matplotlib.pyplot as plt
import os

# E. coli MinE sequence
wtMinE = read_fasta('./Data/EColi_MinE_MSA.fasta')

# sampled seqs
sampled_minEs = read_fasta('../Data/samples.fasta')

# make filter
filter = Filter()

### Filtering step 1: Similarity to E. coli MinE

all_distances, filtered = filter.filter_by_target_identity(wtMinE[list(wtMinE.keys())[0]], sampled_minEs, 0.6)
print('Keeping %i sequences after step 1. ' % (len(filtered)))

# plot distribution of identity
plt.hist(all_distances, bins=20, alpha=0.8, edgecolor = 'grey', color='lightblue')
plt.title('Identitiy to E. coli MinE distribution of sampled sequences')
plt.xlabel('identity')
plt.savefig('../Data/identity_hist.png')

# save new set of sequences to file. From now on this is not a MSA any more!!
filterstep1_fasta_msa = '\n'.join('>' + head + '\n' + seq for head, seq in filtered.items())
filterstep1_path_msa = '../Data/samples_filterstep1_msa.fasta'
open(filterstep1_path_msa, 'w').write(filterstep1_fasta_msa)
filterstep1_fasta = '\n'.join('>' + head + '\n' + seq.replace('-', '') for head, seq in filtered.items())
filterstep1_path = '../Data/samples_filterstep1.fasta'
open(filterstep1_path, 'w').write(filterstep1_fasta)

### Filtering step 1: Clustering by cd-hit
filterstep2_path = '.../Data/samples_filterstep2'
filtered2 = filter.filter_by_cluster_identity(filterstep1_path, filterstep2_path, 0.6, 4)
print('Keeping %i sequences after step 2' % (len(filtered2.keys())))
# msa
filterstep2_fasta_msa = '\n'.join('>' + head + '\n' + seq for head, seq in filtered.items() if head in filtered2.keys())
filterstep2_path_msa = '../Data/samples_filterstep2_msa.fasta'
open(filterstep2_path_msa, 'w').write(filterstep2_fasta_msa)
# each sequence with wildtype MinD for alphafold interaction prediction, in own directory
af_dir_path = '../Data/AlphaFold'
if not os.path.isdir(af_dir_path):
    os.mkdir(af_dir_path)
# read in MinD
minD = read_fasta('./Data/EColi_MinD.fasta')
for head, seq in filtered2.items():
    content = '>' + head + '\n' + seq + '\n>' + list(minD.keys())[0] + '\n' + list(minD.values())[0]
    open(os.path.join(af_dir_path, head + '.fasta'), 'w').write(content)

'''
First: Add M at start of sequence if not there (script 02_addMs.py
Then: Two jobs to do
1. Use AlphaFold to predict interactions using the fasta files stored in ../Data/AlphaFold/
For the next step, put the results in ../Data/AlphaFold_Results/
'''