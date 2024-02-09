import shutil, os, subprocess
from Filtering.Structure import StructurePrediction
from Data.DataUtils.Utils import read_fasta
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import PDBParser
import pandas as pd

'''
Scoring the sampled and clustered sequences by solubility, MinD-interaction, dimerisation and membrane-binding 
capability to find the best synMinE candidates. 
'''


### necessary paths

# path to fasta file with sample sequences
samples_path = '../Data/final_samples.fasta'
# path to root where resulting scores should be saved
results_root = '../Data'
# path to directorys storing AlphaFold PDB files
af_root_minD = '../Data/AlphaFold_MinD_PDBs'
af_root_dimer = '../Data/AlphaFold_Dimer_PDBs'
af_root_npp = '../Data/NPP_Predictions'
# path to directory storing synMinE-MinD interaction PAE matrices
pae_root_minD = '../Data/AlphaFold_MinD_PAEs'
pae_root_dimer = '../Data/AlphaFold_Dimer_PAEs'
pae_root_npp = pae_root_minD # NPP ratios were calculated using PDBs from MinD-interaction
# paths to where PAE images should be saved
pae_image_root_minD = '../Data/AlphaFold_MinD_PAE_Images'
pae_image_root_dimer = '../Data/AlphaFold_Dimer_PAE_Images'
fig_root = '../Data/Hydrophobicity_Traces'



### (1) Solubility
# Using protein-sol to predict the solubilities of the generated and pre-filtered sequences stored in the
# fasta file <samples_path>

# copy file to correct directory
print('Preparing protein-sol...')
shutil.copyfile(samples_path, 'Software/protein-sol-sequence-prediction-software/samples.fasta')

# change working directory to Software/protein-sol-sequence-prediction-software/
# The software has to be downloaded from https://protein-sol.manchester.ac.uk/ !!
os.chdir('Software/protein-sol-sequence-prediction-software/')

# run protein-sol
print('Running protein-sol. This might take a while...')
subprocess.call(['./multiple_prediction_wrapper_export.sh', './samples.fasta'])
print('Done.')

# go to previous wd
os.chdir('../..')

# read in result
print('Parsing result to nice format...')
results = open('Software/protein-sol-sequence-prediction-software/seq_prediction.txt').read().strip().split('\n\n')[1:]
parsed = {}
for result in results:
    name = result.split('\n')[0].split(',')[1][1:]
    value = float(result.split('\n')[0].split(', ')[1])
    parsed[name] = value
print('Values found for %i sequences.' % (len(parsed.keys())))

# save in tsv format
solubilities_path = os.path.join(results_root, 'solubilities.tsv')
tsv = 'Name\tSolubility score\n'
tsv += '\n'.join([name + '\t' + str(value) for name, value in parsed.items()])
open(solubilities_path, 'w').write(tsv)

# make and save histogram
plt.hist(list(parsed.values()), bins=20, alpha=0.8, edgecolor = 'grey', color='lightblue')
plt.title('Predicted solubility distribution')
plt.xlabel('Solubility value')
plt.savefig(os.path.join(results_root, 'solubility_distribution.png'))
plt.clf()

### (2) Interaction with MinD
# Use PAE matrix of AlphaFold dimer prediction of synMinE and MinD (E. coli) to evaluate interaction.
# Also save secondary structure of all synMinEs in "dimer conformation"


# read in sampled sequences
sequences = read_fasta(samples_path)

# make paths to PDB and PAE files
all_pdb_paths = {pdb[:-4].replace('sample', 'Sample_sequence'): os.path.join(af_root_minD, pdb) for pdb in os.listdir(af_root_minD)}
all_pae_paths = {pdb[:-4].replace('sample', 'Sample_sequence'): os.path.join(pae_root_minD, pdb) for pdb in os.listdir(pae_root_minD)}

# storing all values that need to be printed to files later
all_means = {}
sec_structs = {}

# calculate (normalised) mean PAE per alpha helix for MinD-synMinE interaction, get lowest one per sample
print('Calculating scores')
for name, pdb_path in all_pdb_paths.items():
    print(name)

    # prerequisites
    sequence = sequences[name]

    # read in structure
    structure = StructurePrediction(pdb_path, all_pae_paths[name])
    structure.extract_secondary_structure(len(sequence))
    structure.extract_secondary_structure_motifs()

    sec_structs[name] = structure.secondary_structure

    # calculate average PAE to MinD for each alpha helix
    helix_means = []
    for helix in structure.helices:
        helix_mean = np.mean(structure.pae[helix[0]:helix[1], len(sequence):])
        helix_means.append(helix_mean)

    # save lowest result
    all_means[name] = min(helix_means)

# plotting and saving plot
plt.hist(all_means.values(), bins=20, alpha=0.8, edgecolor = 'grey', color='lightblue')
plt.title('Best minD interaction PAE per helix')
plt.xlabel('PAE')
plt.savefig(os.path.join(results_root, 'minD-interaction_paes_hist.png'))
plt.clf()

# saving interaction score data
print('Saving Interaction Data')
interaction_minD_path = os.path.join(results_root, 'interaction_scores.tsv')
interaction_file = 'Name\tPAE Means\n'
interaction_file += '\n'.join([name + '\t' + str(all_means[name]) for name in all_means.keys()])
open(interaction_minD_path, 'w').write(interaction_file)

# saving secondary structure of synMinE
print('Saving Secondary Structure')
sec_str_file = 'Name\tSecondaryStructure\nName\tSequence\n'
sec_str_file += '\n'.join([name + '\t' + str(sec_structs[name]) + '\n' + name + '\t' + str(sequences[name]) for name in sec_structs.keys()])
open(os.path.join(results_root, 'secondary_structures.tsv'), 'w').write(sec_str_file)

### (3) Dimerisation
# Use PAE matrix of AlphaFold homodimer prediction of synMinEs evaluate dimerisation capability.

# make paths to PDB and PAE files
all_pdb_paths = {pdb[:-4].replace('sample', 'Sample_sequence').replace('_dimer', ''): os.path.join(af_root_dimer, pdb) for pdb in os.listdir(af_root_dimer)}
all_pae_paths = {pdb[:-4].replace('sample', 'Sample_sequence').replace('_dimer', ''): os.path.join(pae_root_dimer, pdb) for pdb in os.listdir(pae_root_dimer)}

# storing all values that need to be printed to files later
all_means = {}

# calculate (normalised) mean PAE per alpha helix for MinD-synMinE interaction, get lowest one per sample
print('Calculating scores')
for name, pdb_path in all_pdb_paths.items():
    print(name)

    # prerequisites
    sequence = sequences[name]

    # read in structure
    structure = StructurePrediction(pdb_path, all_pae_paths[name])
    structure.extract_secondary_structure(len(sequence))
    structure.extract_secondary_structure_motifs()

    # calculate mean PAE between different synMinEs per structured region (alpha helix or beta sheet)
    structured_regions = structure.helices + structure.sheets
    structured_regions_means = []
    for region in structured_regions:
        for region_2 in structured_regions:
            region_mean = np.mean(structure.pae[region[0]:region[1], len(sequence) + region_2[0]: len(sequence) + region_2[1]])
            structured_regions_means.append(region_mean)
    final_mean = np.mean(structured_regions_means)

    # calculate and store final values
    all_means[name] = final_mean

# plotting histogram of corr coeffs
plt.hist(all_means.values(), bins=20, alpha=0.8, edgecolor = 'grey', color='lightblue')
plt.title('mean PAEs between structured regions')
plt.xlabel('PAE')
plt.savefig(os.path.join(results_root, 'dimer_PAE_hist.png'))
plt.clf()

# saving corr coeffs
print('Saving PAE Data')
dimerisation_path = os.path.join(results_root, 'dimer_PAEs.tsv')
interaction_file = 'Name\tmean PAEs between structured regions\n'
interaction_file += '\n'.join([name + '\t' + str(all_means[name]) for name in all_means.keys()])
open(dimerisation_path, 'w').write(interaction_file)

### (4) Membrane binding
# Use NPP ration from protein-sol patches (https://protein-sol.manchester.ac.uk/patches, Hebditch M and Warwicker J.Web-based display of protein surface and pH-dependent properties for assessing the developability of biotherapeutics. Scientific Reports (2019))
# at N-terminus to score membrane binding

if not os.path.isdir(fig_root):
    os.mkdir(fig_root)

# make paths to PDB and PAE files
all_pdb_paths = {pdb[:-4].replace('sample', 'Sample_sequence').replace('_ratio', ''): os.path.join(af_root_npp, pdb) for pdb in os.listdir(af_root_npp)}
all_pae_paths = {pdb[:-4].replace('sample', 'Sample_sequence'): os.path.join(pae_root_npp, pdb) for pdb in os.listdir(pae_root_npp)}

# storing all values that need to be printed to files later
all_hydrophobicities = {}
all_n_term_helices = {}

# calculate (normalised) mean PAE per alpha helix for MinD-synMinE interaction, get lowest one per sample
print('Calculating membrane binding scores')
for name, pdb_path in all_pdb_paths.items():
    print(name)

    # prerequisites
    sequence = sequences[name]

    # read in structure
    structure = StructurePrediction(pdb_path, all_pae_paths[name])
    structure.extract_secondary_structure(len(sequence))
    structure.extract_secondary_structure_motifs()

    # parse PDB file
    parser = PDBParser()
    pdb = parser.get_structure(name, structure.pdb_path)

    # calculate hydrophobicity trace for synMinE (average over all atoms per AA)
    hydrophobicity_trace = []
    for i, res in enumerate(pdb[0]['A'].get_residues()):
        hydrophobicity = 0
        n_atoms = 0
        for a in res.get_atoms():
            n_atoms += 1
            hydrophobicity += a.bfactor
        hydrophobicity_trace.append(hydrophobicity/n_atoms)

    # plot trace (TODO nicer layout?)
    plt.plot(hydrophobicity_trace)
    x = np.arange(len(sequence))
    plt.xticks(x, structure.secondary_structure)
    plt.title(name)
    plt.savefig(os.path.join(fig_root, name + '.png'))
    plt.clf()

    # get hydrophobicity score
    hydrophobicity_score = 0
    score_done = False
    if structure.helices != []:
        first_helix = structure.helices[0]
        # check if helix is N-terminal and short
        if first_helix[0] < 10 and first_helix[1] < 30:
            hydrophobicity_score = np.mean(hydrophobicity_trace[first_helix[0]:first_helix[1]]) # if alpha helix: only in helix region
            score_done = True
            all_n_term_helices[name] = True
    if not score_done:
        all_n_term_helices[name] = False
        if structure.sheets != [] and structure.helices != []:
            first_region_start = min(structure.helices[0][0], structure.sheets[0][0])
        elif structure.helices != []:
            first_region_start = structure.helices[0][0]
        elif structure.sheets != []:
            first_region_start = structure.sheets[0][0]
        else:
            hydrophobicity_score = 0
            continue
        hydrophobicity_score = np.mean(hydrophobicity_trace[:first_region_start])

    # calculate and store final values
    all_hydrophobicities[name] = hydrophobicity_score

# plotting histogram of hydrophobicity scores
plt.hist(all_hydrophobicities.values(), bins=20, alpha=0.8, edgecolor = 'grey', color='lightblue')
plt.title('hydrophobiciy scores of N-terminal helix/region')
plt.xlabel('hydrophobicity score')
plt.savefig(os.path.join(results_root, 'membrane_binding_hist.png'))
plt.clf()

# saving corr coeffs
print('Saving Hydrophobicity Data')
membrane_binding_path = os.path.join(results_root, 'membrane_binding_scores.tsv')
interaction_file = 'Name\thydrophobiciy scores of N-terminal helix/region\tn-terminal alpha helix\n'
interaction_file += '\n'.join([name + '\t' + str(all_hydrophobicities[name]) + '\t' + str(all_n_term_helices[name]) for name in all_hydrophobicities.keys()])
open(membrane_binding_path, 'w').write(interaction_file)

### (5) Merge scores
# Normalise all scores to be between 0 and 1, then sum up, sort, and get best 24 and worst 24 seqs.

def normalise(data):
    '''
    normalise a dataset such that data is between 0 and 1
    '''
    return (data - np.min(data)) / (np.max(data) - np.min(data))

# reading in data
solubilities = pd.read_csv(solubilities_path, delimiter='\t')
interaction_minD = pd.read_csv(interaction_minD_path, delimiter='\t')
dimerisation = pd.read_csv(dimerisation_path, delimiter='\t')
membrane_binding = pd.read_csv(membrane_binding_path, delimiter='\t')

# for membrane binding: Add 0.5 if n-terminal region forms proper helix
membrane_binding.loc[membrane_binding['n-terminal alpha helix'] == True, 'hydrophobiciy scores of N-terminal helix/region'] += 0.5

# make all scores between 0 and 1 ("normalise") where 0 is worst and 1 is best
# solubility is in wanted format by default
interaction_minD['PAE Means'] = -(normalise(interaction_minD['PAE Means']) - 1) # large PAE is bad, so reverse direction
dimerisation['mean PAEs between structured regions'] = -(normalise(dimerisation['mean PAEs between structured regions']) - 1) # large PAE is bad, so reverse direction
membrane_binding['hydrophobiciy scores of N-terminal helix/region'] = normalise(membrane_binding['hydrophobiciy scores of N-terminal helix/region'])

# make dfs to dicts
solubilities = pd.Series(solubilities['Solubility score'].values,index=solubilities['Name']).to_dict()
interaction_minD = pd.Series(interaction_minD['PAE Means'].values,index=interaction_minD['Name']).to_dict()
dimerisation = pd.Series(dimerisation['mean PAEs between structured regions'].values,index=dimerisation['Name']).to_dict()
membrane_binding = pd.Series(membrane_binding['hydrophobiciy scores of N-terminal helix/region'].values,index=membrane_binding['Name']).to_dict()

# make final scores
def score(name):
    return solubilities[name] + interaction_minD[name] + dimerisation[name] + membrane_binding[name]

final_scores = {name: score(name) for name in solubilities.keys()}

plt.hist(final_scores.values(), bins=30)
plt.savefig(os.path.join(results_root, 'final_hist.png'))
plt.clf()

open(os.path.join(results_root, 'final_scores.tsv'), 'w').write('\n'.join([name + '\t' + str(value) for name, value in final_scores.items()]))
print('Done with all. Now start experiments.')