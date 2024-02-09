import os, yaml
from Optimiser.Losses import elbo, reconstruction_loss
from torch.utils.data import DataLoader
import torch.optim as optim
import torch
import numpy as np

class Run():

    def __init__(self, config,  train_data, val_data, model, config_root='Configs/Optimiser', silent = False, device = 'cpu', log_nth = 10):

        config_root += '/'
        self.config_root = config_root
        self.config = self.load_config(config)

        self.SILENT = silent
        self.model = model
        self.device = device
        self.log_nth = log_nth

        self.optimiser_config = self.config['optimiser']
        self.loss_func = self.config['loss_func']

        # make dataloader
        self.train_dl = DataLoader(train_data, batch_size=int(self.optimiser_config['batch_size']), shuffle=True)
        self.val_dl = DataLoader(val_data, batch_size=int(self.optimiser_config['batch_size']), shuffle=True)

        # make optimiser
        if self.optimiser_config['optim'] == 'Adam':
            self.optimiser = optim.Adam(self.model.parameters(), lr=self.optimiser_config['learning_rate'])
        else:
            raise NotImplementedError

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

    def run(self):
        '''

        '''

        # if gpu available move to gpu here
        device = torch.device(self.device)
        self.model.to(device)

        # storing history
        train_loss_history = []
        val_loss_history = []

        iterations_per_epoch = len(self.train_dl)
        for epoch in range(self.optimiser_config['num_epochs']):

            if not self.SILENT:
                print('[EPOCH %i] Starting Training...' % (epoch+1))

            for i, (inputs) in enumerate(self.train_dl):
                inputs = inputs.to(device)

                # forward/backward pass
                self.optimiser.zero_grad()
                reconstructions, means, log_vars = self.model(inputs)

                # if specified in configs, only start elbo loss after some lag time
                if epoch <= self.loss_func['cutoff_epoch'] - 1:
                    loss = reconstruction_loss(reconstructions, inputs)
                else:
                    loss = elbo(self.loss_func['error_relation'], inputs, reconstructions, means, log_vars)
                loss.backward()
                self.optimiser.step()

                # logging loss
                train_loss_history.append(loss.cpu().detach().numpy())

                # logging if wanted
                if (not self.SILENT) and self.log_nth and i % self.log_nth == 0 and i != 0:
                    last_nth_train_losses = train_loss_history[-self.log_nth:]
                    cur_mean_loss = np.mean(last_nth_train_losses)
                    print('\r[Iteration %i/%i] TRAIN loss: %.3f' % (
                        i + epoch * iterations_per_epoch,
                        iterations_per_epoch * self.optimiser_config['num_epochs'],
                        cur_mean_loss
                    ), end='')

            # validation. one per epoch
            if not self.SILENT:
                print('\n[EPOCH %i] Training finished, starting Validation...' % (epoch+1))
            val_losses = []
            self.model.eval() # setting model to evaluation mode

            # Check means and stds of latent space
            all_means = 0
            all_deviations = 0
            n_means = 0
            all_stds = 0
            n_stds = 0

            for i, inputs in enumerate(self.val_dl):
                inputs = inputs.to(device)
                reconstructions, means, log_vars = self.model(inputs)

                # for info
                all_means += torch.mean(means)
                all_deviations += torch.mean(torch.abs(means))
                all_stds += torch.mean(torch.exp(0.5*log_vars))
                n_means += 1
                n_stds += 1

                if epoch <= self.loss_func['cutoff_epoch'] - 1:
                    loss = reconstruction_loss(reconstructions, inputs)
                else:
                    loss = elbo(self.loss_func['error_relation'], inputs, reconstructions, means, log_vars)

                loss = loss.cpu().detach().numpy()
                val_losses.append(loss)
            val_loss = np.mean(val_losses)
            val_loss_history.append(val_loss)

            # logging if wanted
            if not self.SILENT:
                print('[EPOCH %i] Validation loss: %.3f, average mean: %.3f, average mean deviation from 0: %.3f,average std: %.3f' % (epoch+1, val_loss, all_means / n_means, all_deviations/ n_means, all_stds / n_stds))

            self.model.train() # setting to training mode again


        self.model.eval()
        return train_loss_history, val_loss_history, self.model