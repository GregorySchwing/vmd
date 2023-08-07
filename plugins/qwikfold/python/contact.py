import os 
import sys
from matplotlib import pyplot as plt
import jax
import pickle

filename=sys.argv[1]

with open(filename,'rb') as f :
    prediction_result = pickle.load(f)  

dist_bins   = jax.numpy.append(0,prediction_result["distogram"]["bin_edges"])
dist_mtx    = dist_bins[prediction_result["distogram"]["logits"].argmax(-1)]

dist=jax.nn.softmax(prediction_result["distogram"]["logits"])[:,:,dist_bins < 8].sum(-1)
plt.title(f"Contact Map\n{os.path.basename(filename)}")

plt.imshow(dist,cmap='binary')
plt.show()
