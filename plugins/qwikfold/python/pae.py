import sys
import os
from matplotlib import pyplot as plt
import jax
import pickle

filename=sys.argv[1]

with open(filename,'rb') as f :
    prediction_result = pickle.load(f)  

dist=prediction_result["predicted_aligned_error"]
plt.title(f"Predicted Alignment Error (PAE)\n{os.path.basename(filename)}")
plt.imshow(dist,cmap='bwr')
plt.colorbar()
plt.show()

