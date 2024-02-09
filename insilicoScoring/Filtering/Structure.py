import subprocess
import shutil
import numpy as np

# stride MUST be installed!
assert shutil.which('stride'), "Looks like stride is not installed. Stride MUST be installed for this analysis to run. Install from http://webclu.bio.wzw.tum.de/stride/install.html"

class Stride():
    def __init__(self, stride_string):
        self.string = stride_string

    def extract_secondary_structure(self, string):
        return ''.join([x[10:60] for x in string.split('\n') if x[:3] == 'STR'])




class StructurePrediction():
    def __init__(self, pdb_path, pae_path):
        self.pdb_path = pdb_path
        self.pae = self.read_pae(pae_path)
        self.secondary_structure = None

    def extract_secondary_structure(self, end=None):
        p = subprocess.run(['stride', '-o', self.pdb_path], stdout=subprocess.PIPE)
        stride_output = p.stdout.decode('utf-8')
        s = Stride(stride_output)
        sec_str = s.extract_secondary_structure(s.string)
        if end == None:
            self.secondary_structure = sec_str
        else:
            self.secondary_structure = sec_str[:end]

    def read_pae(self, pae_path):
       return np.genfromtxt(pae_path, delimiter=',')

    def extract_secondary_structure_motifs(self):
        # find alpha helices
        helices = []
        sheets = []
        is_helix = False
        is_sheet = False
        for i in range(len(self.secondary_structure)):
            if is_helix == False and self.secondary_structure[i] == 'H':
                is_helix = True
                start = i
            elif is_helix == True and self.secondary_structure[i] != 'H':
                is_helix = False
                end = i
                helices.append([start, end])
            elif is_sheet == False and self.secondary_structure[i] == 'E':
                is_sheet = True
                start = i
            elif is_sheet == True and self.secondary_structure[i] != 'E':
                is_sheet = False
                end = i
                sheets.append([start, end])
        self.helices = helices
        self.sheets = sheets