import sys
import os
from matplotlib import pyplot as plt
import jax
import pickle

filename=sys.argv[1]

with open(filename,'rb') as f :
    prediction_result = pickle.load(f)  

dist=prediction_result["distogram"]["logits"].argmax(-1)
plt.title(f"Distance Map\n{os.path.basename(filename)}")
plt.imshow(dist,cmap='viridis_r')
plt.colorbar()
plt.show()

