import numpy as np
import torch

def read_fasta(path):
    '''
    Reading in a fasta file
    :param path: path to fasta file to read in
    :return: dictionary {header: seq}
    '''
    entries = open(path).read().strip().split('>')[1:]
    return {entry.split('\n')[0]: ''.join(entry.split('\n')[1:]) for entry in entries}

def one_hot_encode(sequence):
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
    encoding = np.zeros(( 21*len(sequence)))
    for i, aa in enumerate(sequence):
        if not aa in AAs:
            print('Unknown amino acid %s in protein. Not encoded. ' % (aa))
            return -1
        encoding[ i*21 + AAs[aa]] = 1
    return encoding

def make_sequence(one_hot_encoding):
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