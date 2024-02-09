from random import shuffle
import random
from Data.DataUtils.Utils import read_fasta
'''
rankes the proteins by final scores, gives proteins random IDs and saves them to fasta file
'''

# read in data
sequences = read_fasta('../Data/final_samples.fasta')
scores_path = '../Data/final_scores.tsv'
proteins = {line.split('\t')[0]: float(line.split('\t')[1]) for line in open(scores_path).read().strip().split('\n')}

# sort by score
proteins = dict(sorted(proteins.items(), key=lambda item: item[1]))

# select best 24
best_24_names = list(proteins.keys())[-24:]

# select worst 24
worst_24_names = list(proteins.keys())[:24]

# mix randomly and assign names
random.seed(42)
mixed_selected = best_24_names + worst_24_names
shuffle(mixed_selected)
map = {name: 'synMinEv' + str(i+1) for i, name in enumerate(mixed_selected)}

# save map as file
open('../Data/blinding_map.tsv', 'w').write('\n'.join([name + '\t' + new_name for name, new_name in map.items()]))

# save fasta file with all sequences with new names
new_fasta = {new_name: sequences[name] for name, new_name in map.items()}
open('../Data/synMinEs.fasta', 'w').write('\n'.join(['>' + new_name + '\n' + seq for new_name, seq in new_fasta.items()]))

print('Done.')
