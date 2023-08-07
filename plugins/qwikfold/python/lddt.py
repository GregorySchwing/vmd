import os
import sys
from matplotlib import pyplot as plt
import pickle

filename=sys.argv[1]

with open(filename,'rb') as f :
    prediction_result = pickle.load(f)  

data=prediction_result['plddt']

plt.plot(data)
plt.title(f" Local Distance Difference Test (lDDT)\n{os.path.basename(filename)}")
plt.ylabel('lDDT\n(predicted)')
plt.xlabel('Residue')
plt.show()
