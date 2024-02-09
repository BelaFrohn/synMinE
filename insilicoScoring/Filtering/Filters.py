from Data.DataUtils.Clustering import cdhit
from Data.DataUtils.Utils import read_fasta

class Filter():
    def __init__(self):
        pass

    def filter_by_target_identity(self, target_seq, seqs, threshold):
        '''
        Filters out sequences that are identical in too many positions.
        Identity is calculated by identical aas / min(length(A),length(B))
        Sequences are kept if identity < threshold.
        :param target_seq: sequence to compare to
        :param seqs: sequences to filter, as dict
        :param threshold: float
        :return: filtered dict
        '''
        outset = {}
        all_distances = []
        length_target = len(target_seq.replace('-', ''))
        for head, seq in seqs.items():
            length_current = len(seq.replace('-', ''))
            identical = 0
            for i, aa in enumerate(seq):
                if aa == '-' or target_seq[i] == '-':
                    continue
                if aa == target_seq[i]:
                    identical += 1
            identity = identical / min(length_target, length_current)
            all_distances.append(identity)
            if identity < threshold:
                outset[head] = seq
        return all_distances, outset

    def filter_by_cluster_identity(self, in_path, out_path, threshold, wordlength=5):
        '''
        uses cd-hit to cluster sequences by identity.
        :param seqs: sequences to filter, as dict
        :param threshold: float
        :return: filtered dict
        '''
        clusterer = cdhit()
        clusterer.c = threshold
        clusterer.n = wordlength
        clusterer.cluster(in_path, out_path)
        return read_fasta(out_path + '.fasta')