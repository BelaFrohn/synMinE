import torch
import torch.nn as nn

def reconstruction_loss(reconstruction, original):
    loss = nn.BCELoss()
    l = loss(reconstruction, original)
    return l

def kl_divergence(mean, log_var): # kullback leibler divergence (loss func for distribution)
    loss = torch.sum(-0.5 * (1 + log_var - mean ** 2 - log_var.exp())) / (mean.size()[0] * mean.size()[1] * mean.size()[2])
    return loss

def elbo(error_relation, original, reconstruction, means, log_vars):
    rec = reconstruction_loss(reconstruction, original)
    kl = kl_divergence(means, log_vars)
    elbo = error_relation*kl + rec
    return torch.mean(elbo)