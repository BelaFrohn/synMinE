from Data.dataset import Dataset
from Models.FC import fcVAE
from Optimiser.optimiser import Run
import matplotlib.pyplot as plt
from Distribution.Utils import sample2, project_to_latent
import numpy as np
import json
import torch

# define configs
data_config = 'chosen'
model_config = 'chosen'
optimiser_config = 'chosen'

# all rest is automatised

# loading data
print('Loading training data...')
train_data = Dataset(data_config, 'Train')
print('Training data with', len(train_data), 'datapoints.')
print('Loading training data...')
val_data = Dataset(data_config, 'Val')
print('Validation data with', len(val_data), 'datapoints.')

# create model
print('Creating model.')
model = fcVAE(model_config)

# count parameters
param_count = 0
for p in model.parameters():
    param_count += p.numel()
print('Model with', param_count, 'parameters.')

# make run
run = Run(optimiser_config, train_data, val_data, model)
result = run.run()

# print training curve
x1 = np.linspace(1,len(train_data) * run.optimiser_config['num_epochs'], len(result[0]))
x2 = np.linspace(1,len(train_data) * run.optimiser_config['num_epochs'], len(result[1]))
plt.plot(x1, result[0])
plt.plot(x2, result[1])
plt.savefig('./Results_rerun/vae_training.png')

# sample from latent space
print('Random samples from latent space')
true_path = train_data.paths['filter_postMSA'] + '.fasta'
print(true_path)
sample_fasta, frq_corr, cofrq_corr = sample2(result[2], true_path, model.hparams['latent_dims'] , 4000)
print('AA frequency correlation:', frq_corr)
print('AA pairwise frequency correlation:', cofrq_corr)
open('./Results_rerun/samples.fasta', 'w').write(sample_fasta)
print('Saved sample fasta file.')

# save training data to make nicer plot later
#open('./Results_rerun/vae_losses.json', 'w').write(json.dumps(result.tolist()))

# save model
torch.save(model.state_dict(), './Results_rerun/model')

# get and save latent values of train data
projections = {}

print('Reading in full dataset.')
full_data = Dataset(data_config, 'Full')

print('Calculating projections.')
for i, x in enumerate(full_data):
    cur_projections = project_to_latent(result[2], x)
    print(cur_projections)
    id = full_data.data_ids[i]
    projections[id] = cur_projections

print('Number of projections: %i' % (len(projections)))

projections_path = './Results_rerun/projections.txt'

print('Writing projections to file')
txt = '\n'.join([id + '\t' + str(projection) for id, projection in projections.items()])
open(projections_path, 'w').write(txt)





