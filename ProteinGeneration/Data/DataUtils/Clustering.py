import subprocess
import os
from Data.DataUtils import Utils

class cdhit:
    def __init__(self):
        # identity threshold
        self.c = 0.9
        # word length
        self.n = 5
        # memory limit
        self.M = 32000
        # length of description in .clstr file
        self.d = 200
        # number of threads (0 = all)
        self.T = 32

    def cluster(self, input, output):
        # call
        command = 'cd-hit -i ' + input + ' -o ' + output + ' -c ' + str(self.c) + ' -n ' + str(self.n) + ' -M ' + str(int(self.M)) + ' -d ' + str(self.d) + ' -T ' + str(int(self.T))
        subprocess.call(command, shell=True)
        # make fasta file to fasta file
        os.rename(output, output + '.fasta')

class IdentitySearch():
    '''
    My own version of clustering, simply clusters by 100% identity
    '''
    def __init__(self):
        self.SILENT = False

    def cluster(self, input, output):
        in_seqs = Utils.read_fasta(input)
        clusters = {}
        for header, seq in in_seqs.items():
            if seq in clusters.keys():
                clusters[seq].append(header)
            else:
                clusters[seq] = [header]
        if not self.SILENT:
            print('Grouped into', len(clusters.keys()), 'identical clusters.')
        # write fasta file
        fasta_content = '\n'.join(['>' + clusters[seq][0] + '\n' + seq for i, seq in enumerate(list(clusters.keys()))])
        open(output + '.fasta', 'w').write(fasta_content)
        # write cluster file
        cluster_content = '\n'.join(['>cluster_' + str(i) + '\n' + '\n'.join([str(j) + '\t' + str(len(seq)) + 'aa, >' + header for j, header in enumerate(clusters[seq])]) for i, seq in enumerate(list(clusters.keys()))])
        open(output + '.clstr', 'w').write(cluster_content)
        if not self.SILENT:
            print('Saved files.')




