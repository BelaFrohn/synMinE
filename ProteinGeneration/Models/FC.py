import torch
import torch.nn as nn
import os
import yaml
import numpy as np


class fcVAE(nn.Module):
    def __init__(self, config, config_root='Configs/Model'):
        super().__init__()

        config_root += '/'
        self.config_root = config_root

        self.config = self.load_config(config)

        self.architecture = self.config['architecture']
        self.hparams = self.config['hparams']

        # building encoder
        encoding_layers = []
        for i in range(self.architecture['layers']):
            if i == 0:
                encoding_layers += self.make_layer(self.architecture['input_size'], self.hparams['hidden_size_1'], self.architecture['activation_function'])
            else:
                encoding_layers += self.make_layer(self.hparams['hidden_size_' + str(i)], self.hparams['hidden_size_' + str(i+1)], self.architecture['activation_function'])
            if self.architecture['dropout']:
                encoding_layers += [nn.Dropout()]
            if self.architecture['batchnorm']:
                encoding_layers += [nn.BatchNorm1d(1)]
        self.encoder = nn.Sequential(*encoding_layers)

        # building middle part
        self.mean = nn.Sequential(
            nn.Linear(self.hparams['hidden_size_' + str(self.architecture['layers'])], self.hparams['latent_dims']),
        )
        self.log_variance = nn.Sequential(
            nn.Linear(self.hparams['hidden_size_' + str(self.architecture['layers'])], self.hparams['latent_dims']),
        )

        # building decoder
        decoding_layers = []
        for i in range(self.architecture['layers'], 0, -1):
            if i == self.architecture['layers']:
                decoding_layers += self.make_layer(self.hparams['latent_dims'], self.hparams['hidden_size_' + str(i)], self.architecture['activation_function'])
            else:
                decoding_layers += self.make_layer(self.hparams['hidden_size_' + str(i+1)], self.hparams['hidden_size_' + str(i)], self.architecture['activation_function'])
            if self.architecture['dropout']:
                decoding_layers += [nn.Dropout()]
            if self.architecture['batchnorm']:
                decoding_layers += [nn.BatchNorm1d(1)]

        decoding_layers += [nn.Linear(self.hparams['hidden_size_1'], self.architecture['input_size'])]

        self.decoder = nn.Sequential(*decoding_layers)


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

    def reparameterize(self, means, stds):
        rand = torch.randn_like(means)
        return stds * rand + means

    def softmax_on_bins(self, one_hot_encoded):
        stacked = torch.stack(torch.split(one_hot_encoded, 21, 2), 3)
        softmaxed = nn.functional.softmax(stacked, 2)
        flattened = torch.flatten(softmaxed, 2)
        return flattened

    def forward(self, x):

        # encode from one hot encoding to latent space
        y = self.encoder(x)
        means = self.mean(y)
        log_vars = self.log_variance(y)
        stds = torch.exp(0.5*log_vars)

        # reparameterization trick
        latents = self.reparameterize(means, stds)

        # decode from latent space to one hot encoding
        z = self.decoder(latents)

        # softmax per AA
        result = self.softmax_on_bins(z)

        return result, means, log_vars

    def make_layer(self, input_size, output_size, activation_function):
        if activation_function == 'ReLU':
            return [nn.Linear(input_size, output_size), nn.ReLU()]
        else:
            raise NotImplementedError
