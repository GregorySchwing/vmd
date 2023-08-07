import sys
import pathlib
from matplotlib import pyplot as plt
from alphafold.data import parsers
import numpy as np

filename=sys.argv[1]

with open(filename,'r') as f :
    sequences=f.read()
#msa, deletion_matrix, target_names = parsers.parse_stockholm(mgnify)


if pathlib.Path(filename).suffix == '.sto' : 
    MSA = parsers.parse_stockholm(sequences)

if pathlib.Path(filename).suffix == '.a3m' : 
    MSA = parsers.parse_a3m(sequences)


msa = MSA.sequences

aa_map = {restype: i for i, restype in enumerate('ABCDEFGHIJKLMNOPQRSTUVWXYZ-')}
msa_arr = np.array([[aa_map[aa] for aa in seq] for seq in msa])
num_alignments, num_res = msa_arr.shape

fig = plt.figure(figsize=(6, 3))
plt.title('Non-Gap Amino Acids\nin Multiple Sequence Alignment',fontsize=10,y=1)
plt.plot(np.sum(msa_arr != aa_map['-'], axis=0), color='black')
plt.xlabel('Residue')
plt.ylabel('Count')
plt.yticks(range(0, num_alignments + 1, max(1, int(num_alignments / 3))))

fig.tight_layout()

plt.show()
