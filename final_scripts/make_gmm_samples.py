
import matplotlib
matplotlib.use('Agg')
from matplotlib.pyplot import *
import localgroup
import triangle
import sklearn
from sklearn import mixture
import numpy as np
import pickle
import matplotlib.patches as mpatches
import sys

save_path = sys.argv[1]
# Inside the Likelihood object is a "triplet" object called T, which contains an array of sample local groups, each with kinematic parameters consistent with the observational data. Let's plot these kinematic parameters in a "triangle" figure, to show all their 1 and 2-D marginal distributions.

# In[2]:
with open('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/L_Q_presym.pickle', 'rb') as Lfile:
    L = pickle.load(Lfile)

with open('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr_sym_N30_v80_presym.pickle', 'rb') as Trfile:
    Tr = pickle.load(Trfile)

Tr.GMM_sample(200000000, L, reps=1)
with open(save_path, 'wb') as f:
    pickle.dump(Tr.gmm_samples, f)
