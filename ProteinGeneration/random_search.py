from Data.dataset import Dataset
from Models.FC import fcVAE
from Optimiser.optimiser import Run
import matplotlib.pyplot as plt
from Distribution.Utils import sample
import numpy as np
import yaml, json, hashlib, os
import torch

class RandomSearch():
    def __init__(self, name, config, n, n_sample, config_root='Configs/Searches'):

        config_root += '/'
        self.name = name
        self.config_root = config_root
        self.config = self.load_config(config)

        self.n = n
        self.n_sample = n_sample

        self.train_data = Dataset(self.config['data'], 'Train')
        self.val_data = Dataset(self.config['data'], 'Val')



    def load_config(self, config):
        '''
        From TC
        Parse a config name, path, or object to a config object.
        '''
        # config is not a config object
        if type(config) == str:
            # config is a config name
            if not '/' in config and os.path.isfile(self.config_root+config+'.yml'):
                config = self.config_root+config+'.yml'
            # config is a config path
            if os.path.isfile(config):
                with open(config,'r') as file:
                    config = yaml.safe_load(file)
        return config

    def hash(self, config):
        '''
        compute the hash value of a dictionary-type config
        '''
        string = str(config)
        return hashlib.md5(string.encode('utf-8')).hexdigest()

    def run(self, path, model_config, optimiser_config, n_samples):
        '''
        A single training run
        '''

        # check if previously done
        if os.path.exists(path):
            print('This hparam set was already tested')
            return -1

        # make new directory
        os.mkdir(path)

        # make model
        model = fcVAE(model_config)

        # count parameter
        param_count = 0
        for p in model.parameters():
            param_count += p.numel()

        # train
        run = Run(optimiser_config, self.train_data, self.val_data, model)
        run.SILENT = True
        result = run.run()

        # print training curve
        x1 = np.linspace(1,len(self.train_data) * run.optimiser_config['num_epochs'], len(result[0]))
        x2 = np.linspace(1,len(self.train_data) * run.optimiser_config['num_epochs'], len(result[1]))
        plt.plot(x1, result[0])
        plt.plot(x2, result[1])
        plt.savefig(os.path.join(path, 'training.png'))
        plt.clf()

        # calculate statistics to compare
        true_path = self.train_data.paths['filter_postMSA'] + '.fasta'
        sample_fasta, frq, cofrq = sample(result[2], true_path, model.hparams['latent_dims'], n_samples)

        # save sample and statistics
        open(os.path.join(path, 'samples.fasta'), 'w').write(sample_fasta)
        open(os.path.join(path, 'statistics.txt'), 'w').write('AA frequency correlation\t' + str(frq) + '\nAA pairwise frequency correlation\t' + str(cofrq) + '\nParameter count\t' + str(param_count))

        # save model
        torch.save(model, os.path.join(path, 'model.torch'))

        return frq, cofrq

    def search(self):
        '''
        perform random search with n runs.
        '''

        # make n configs
        configs = {}
        while len(configs) < self.n:
            model_config = self.config['model']
            optimiser_config = self.config['optimiser']
            cur_model_architecture = {}
            cur_model_hparams = {}
            cur_optimiser_optim = {}
            cur_optimiser_lossfunc = {}

            for k,v in model_config['architecture'].items():
                if isinstance(v, list):
                    cur_model_architecture[k] = np.random.choice(v)
                else:
                    cur_model_architecture[k] = v
            for k,v in model_config['hparams'].items():
                if isinstance(v, list):
                    cur_model_hparams[k] = np.random.choice(v)
                else:
                    cur_model_hparams[k] = v

            for k,v in optimiser_config['optimiser'].items():
                if isinstance(v, list):
                    cur_optimiser_optim[k] = np.random.choice(v)
                else:
                    cur_optimiser_optim[k] = v
            for k,v in optimiser_config['loss_func'].items():
                if isinstance(v, list):
                    cur_optimiser_lossfunc[k] = np.random.choice(v)
                else:
                    cur_optimiser_lossfunc[k] = v

            cur_model_config = {'architecture': cur_model_architecture, 'hparams': cur_model_hparams}
            cur_optimiser_config = {'optimiser': cur_optimiser_optim, 'loss_func': cur_optimiser_lossfunc}
            cur_total_config = {'model': cur_model_config, 'optimiser': cur_optimiser_config}
            hash = self.hash(cur_total_config)
            if not hash in configs.keys():
                configs[hash] = cur_total_config

        results = {}
        for hash, config in configs.items():
            print('Running', hash)
            result = self.run('Data/Runs/' + hash, config['model'], config['optimiser'], self.n_sample)
            if result != -1: # -1 if already calculated
                results[hash] = result
                open('Data/Runs/' + hash + '/hparams.json', 'w').write(str(config))

        open('Data/Runs/' + self.name + '_results.txt', 'w').write(json.dumps(results))





