import subprocess
import shutil
'''
import os
from Utils.Utils import load_fasta, make_fasta
from Plotting import HelixPlotter, SheetPlotter
import matplotlib.pyplot as plt
import biotite.sequence as seq
import biotite.sequence.graphics as graphics
'''

class Stride():
    def __init__(self, stride_string):
        self.string = stride_string

    def extract_secondary_structure(self, string):
        return ''.join([x[10:60] for x in string.split('\n') if x[:3] == 'STR'])



# stride MUST be installed!
assert shutil.which('stride'), "Looks like stride is not installed. Stride MUST be installed for this analysis to run. Install from http://webclu.bio.wzw.tum.de/stride/install.html"

def extract_secondary_structure(pdb_path):
    p = subprocess.run(['stride', '-o', pdb_path], stdout=subprocess.PIPE)
    stride_output = p.stdout.decode('utf-8')
    s = Stride(stride_output)
    sec_str = s.extract_secondary_structure(s.string)
    return sec_str

'''

# sec structure analysis
sec_structures = {}
parsed_sec_structures = {}
for i1, prots in enumerate(chmp_classes):
    class_name = chmp_class_names[i1]
    print('Analysing', class_name)
    sec_structures[class_name] = []
    aligned_parsed_sec_structures[class_name] = {}
    maps_to_msa_inds[class_name] = []
    parsed_sec_structures[class_name] = []
    aligned_seqs[class_name] = []
    for i2, prot in enumerate(prots):
        struct_path = chmp_paths[i1] + prot + '/' + prot + '/ranked_0.pdb'
        p = subprocess.run(['stride', '-o', struct_path], stdout=subprocess.PIPE)
        stride_output = p.stdout.decode('utf-8')
        s = Stride(stride_output)
        sec_str = s.extract_secondary_structure(s.string)
        sec_structures[class_name].append(sec_str)
        parsed_sec_str = sec_str.replace('T', 'X').replace('G', 'X').replace('B', 'X').replace(' ', 'X')
        aligned_seq = msa[prot]
        aligned_parsed_sec_str = ''
        i3 = 0
        map_to_msa_indices = {}
        for i4, aa in enumerate(aligned_seq):
            if aa == '-':
                aligned_parsed_sec_str += '-'
            else:
                aligned_parsed_sec_str += parsed_sec_str[i3]
                map_to_msa_indices[i3] = i4
                i3+=1
        maps_to_msa_inds[class_name].append(map_to_msa_indices)
        aligned_parsed_sec_structures[class_name][prot] = aligned_parsed_sec_str
        parsed_sec_structures[class_name].append(parsed_sec_str)
        aligned_seqs[class_name].append(aligned_seq)

    print('Saving aligned secondary structures as fasta, where H means helix and X means no helix. ')
    open('Data/Results/' + class_name + '_aligned_secondary_structures.fasta', 'w').write(make_fasta(aligned_parsed_sec_structures[class_name]))

# plotting
for i1, prots in enumerate(chmp_classes):
    class_name = chmp_class_names[i1]
    print('Plotting', class_name)
    fig = plt.figure(figsize=(8.0, len(prots)*0.8))
    for i2, prot in enumerate(prots):
        is_helix = False
        features = []
        start = 0
        end = 0
        for i, aa in enumerate(parsed_sec_structures[class_name][i2]):
            if is_helix == False and aa == 'H':
                is_helix = True
                start = i
            elif is_helix == True and aa != 'H':
                is_helix = False
                end = i
                if maps_to_msa_inds[class_name][i2][end] < 90:
                    features.append(seq.Feature("SecStr", [seq.Location(start, end)], {"sec_str_type" : "helix", "color": "red"}))
                else:
                    features.append(seq.Feature("SecStr", [seq.Location(start, end)], {"sec_str_type" : "helix", "color": "black"}))

        annotation = seq.Annotation(features)


        ax = fig.add_subplot(len(prots),1,i2+1)
        ax.title.set_text(prot)
        graphics.plot_feature_map(
            ax, annotation, multi_line=False, loc_range=(1,len(aligned_seqs[class_name][i2].replace('-',''))),
            feature_plotters=[HelixPlotter(), SheetPlotter()]
        )
    fig.tight_layout()
    plt.savefig('Data/Results/' + class_name + '.png')
'''