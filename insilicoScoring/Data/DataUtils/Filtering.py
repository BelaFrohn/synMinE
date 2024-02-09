import numpy as np

class SeqFilter():
    def __init__(self):
        self.SILENT = False

    def filter_max_length(self, seqs, max_length):
        '''
        Filters sequences by maximum length. input must be dict {header: seq}, returns same
        '''
        return {head: seq for head, seq in seqs.items() if len(seq) <= max_length}

    def filter_min_length(self, seqs, min_length):
        '''
        Filters sequences by minimum length. input must be dict {header: seq}, returns same
        '''
        return {head: seq for head, seq in seqs.items() if len(seq) >= min_length}

    def filter_standard_aas(self, seqs):
        '''
        Filters out sequences with non-standard AAs. input must be dict {header: seq}, returns same
        '''
        standard_aas = 'ARNDCEQGHILKMFPSTWYV'
        out = {}
        for head, seq in seqs.items():
            ok = True
            for aa in seq:
                if not aa in standard_aas:
                    ok = False
                    break
            if ok:
                out[head] = seq
        return out

    def delete_empty_msa_columns(self, msa):
        '''
        Deletes columns that only contain gaps from a MSA. input must be dict {header: seq}, seqs equal length
        '''
        gap_positions = []
        for i in range(len(list(msa.values())[0])):
            all_gaps = True
            for head, seq in msa.items():
                if seq[i] != '-':
                    all_gaps = False
                    break
            if all_gaps:
                gap_positions.append(i)
        out = {}
        for head, seq in msa.items():
            for i in reversed(gap_positions):
                seq = seq[:i] + seq[i+1:]
            out[head] = seq
        return out

    def delete_seqs_with_rare_positions(self, msa, threshold=0.05):
        '''
        Deletes all sequences that have AAs in columns where less than threshold proteins have AAs.
        Thus reduces width of MSA. msa input must be dict {header: seq}
        '''
        frequencies = np.zeros((len(list(msa.values())[0])))
        for head, seq in msa.items():
            for i, aa in enumerate(seq):
                if aa != '-':
                    frequencies[i] += 1
        frequencies = frequencies / len(msa.values())
        indices = np.where(frequencies < threshold)[0]

        out = {}
        for head, seq in msa.items():
            to_delete = False
            for i in indices:
                if seq[i] != '-':
                    to_delete = True
                    break
            if not to_delete:
                out[head] = seq
        return out







