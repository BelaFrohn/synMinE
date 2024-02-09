import numpy as np
import torch

AAs = {
    'A': 0,
    'C': 1,
    'D': 2,
    'E': 3,
    'F': 4,
    'G': 5,
    'H': 6,
    'I': 7,
    'K': 8,
    'L': 9,
    'M': 10,
    'N': 11,
    'P': 12,
    'Q': 13,
    'R': 14,
    'S': 15,
    'T': 16,
    'V': 17,
    'W': 18,
    'Y': 19,
    '-': 20
}

def calculate_aa_frequencies(seqs):
    '''
    Calculates amino acid frequencies for every position given a set of (aligned) sequences
    :param seqs: list of sequences. All must be of same length.
    :return: frequency matrix
    '''
    frq_matrix = np.zeros((21, len(seqs[0])))
    for seq in seqs:
        for i, aa in enumerate(seq):
            if aa == 'X':
                continue
            frq_matrix[AAs[aa], i] += 1
    frq_matrix = frq_matrix/frq_matrix.sum(axis=0, keepdims=1)
    return frq_matrix

def calculate_pairwise_frqs(seqs):
    '''
    Calculates frequencies of pairs of AAs occuring together
    :param seqs: list of sequences. All must be of same length.
    :return: frequency matrix
    '''
    frq_matrix = np.zeros((21*len(seqs[0]), 21*len(seqs[0])))
    for seq in seqs:
        for i, aa in enumerate(seq):
            if aa == 'X':
                continue
            for j, aa2 in enumerate(seq):
                if aa2 == 'X':
                    continue
                frq_matrix[21*i + AAs[aa], 21*j + AAs[aa2]] += 1
    for i in range(21*len(seqs[0])):
        if frq_matrix[i,i] == 0:
            continue
        frq_matrix[i,:] = frq_matrix[i,:] / frq_matrix[i,i]

    return frq_matrix

def calculate_pairwise_frqs_2(seqs):
    '''
    Calculates frequencies of pairs of AAs occuring together
    :param seqs: list of sequences. All must be of same length.
    :return: frequency matrix
    '''
    frq_matrix = np.zeros((20*len(seqs[0]), 20*len(seqs[0])))
    for seq in seqs:
        for i, aa in enumerate(seq):
            if aa == '-':
                continue
            for j, aa2 in enumerate(seq):
                if aa2 == '-':
                    continue
                frq_matrix[20*i + AAs[aa], 20*j + AAs[aa2]] += 1
    return frq_matrix



def make_sequence(one_hot_encoding):
    '''
    make protein sequence from one hot encoding (encoding with gaps)
    '''
    AAs = {
        0: 'A',
        1: 'C',
        2: 'D',
        3: 'E',
        4: 'F',
        5: 'G',
        6: 'H',
        7: 'I',
        8: 'K',
        9: 'L',
        10: 'M',
        11: 'N',
        12: 'P',
        13: 'Q',
        14: 'R',
        15: 'S',
        16: 'T',
        17: 'V',
        18: 'W',
        19: 'Y',
        20: '-'
    }
    ohe = torch.flatten(one_hot_encoding)
    aa = ''.join([AAs[np.argmax(ohe[i: i+21]).item()] for i in range(0, len(ohe), 21)])
    return aa

def sample(model, true_path, latent_dims, n=1000):
    '''

    :param model: trained model
    :param true_path: path to fasta file with true msa (to compare frequencies)
    :return: msa with generated seqs in fasta format, frq correlation, aa pairwise frq correlation
    '''
    model.eval()


    # true aa frequency table
    true_fasta = [''.join(entry.split('\n')[1:]) for entry in open(true_path).read().strip().split('>')[1:]]
    aa_freq_table_true = calculate_aa_frequencies(true_fasta)
    pw_freq_table_true = calculate_pairwise_frqs(true_fasta)

    with torch.no_grad():

        # sample from latent distribution, 100 samples
        latent = torch.randn(n, 1, latent_dims, device='cpu')

        # reconstruct sequences
        sample_seqs = model.decoder(latent)
        sample_seqs = model.softmax_on_bins(sample_seqs)

        # make fasta from reconstructions
        sample_fasta = ''
        seqs = []
        for i in range(n):
            seq = str(make_sequence(sample_seqs.cpu().detach()[i,:,:]))
            sample_fasta += '>Sample_sequence_' + str(i) + '\n' + seq + '\n'
            seqs.append(seq)

        # aa frequencies
        aa_freq_table = calculate_aa_frequencies(seqs)
        pw_freq_table = calculate_pairwise_frqs(seqs)
        frq_corr = np.corrcoef(aa_freq_table_true.flatten(), aa_freq_table.flatten())[0,1]
        cofrq_corrs = np.corrcoef(pw_freq_table_true.flatten(), pw_freq_table.flatten())[0,1]

        return sample_fasta, frq_corr, cofrq_corrs



def sample2(model, true_path, latent_dims, n=1000):
    '''

    :param model: trained model
    :param true_path: path to fasta file with true msa (to compare frequencies)
    :return: msa with generated seqs in fasta format, frq correlation, aa pairwise frq correlation
    '''
    model.eval()


    # true aa frequency table
    true_fasta = [''.join(entry.split('\n')[1:]) for entry in open(true_path).read().strip().split('>')[1:]]
    n_true = len(true_fasta)
    print('True sample size:', n_true)
    print('Calculate true frq')
    true_frq = calculate_pairwise_frqs_2(true_fasta) / n_true
    print('Calculate true triangle')
    true_triangle = list(true_frq[np.triu_indices(len(20*true_fasta[0]), k=1)])
    aa_freq_table_true = calculate_aa_frequencies(true_fasta)

    with torch.no_grad():

        # sample from latent distribution, n samples
        latent = torch.randn(n, 1, latent_dims, device='cpu')

        # reconstruct sequences
        sample_seas = model.decoder(latent)
        sample_seqs = model.softmax_on_bins(sample_seqs)

        # make fasta from reconstructions
        sample_fasta = ''
        seqs = []
        for i in range(n):
            seq = str(make_sequence(sample_seqs.cpu().detach()[i,:,:]))
            sample_fasta += '>Sample_sequence_' + str(i) + '\n' + seq + '\n'
            seqs.append(seq)

        # aa frequencies
        aa_freq_table = calculate_aa_frequencies(seqs)

        # correlations
        frq_corr = np.corrcoef(aa_freq_table_true.flatten(), aa_freq_table.flatten())[0,1]

        n_sample = len(seqs)
        print('Calculate sample frq')
        sample_frq = calculate_pairwise_frqs_2(seqs) / n_sample

        print('Calculate sample triangle')
        sample_triangle = list(sample_frq[np.triu_indices(len(20*seqs[0]), k=1)])

        print('Calculate correlation')
        cofrq_corr = np.corrcoef(true_triangle, sample_triangle)[0,1]

        return sample_fasta, frq_corr, cofrq_corr


def project_to_latent(model, ohe):
    # projects sequences onto latent dimensions

    model.eval()

    with torch.no_grad():
        y = model.encoder(ohe)
        means = model.mean(y)

    return means

















