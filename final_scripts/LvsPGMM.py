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

save_path = "/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/final_scripts/"
Lfile = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/L.pickle'
Trfile = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr.pickle'

with open(Lfile, 'rb') as Lf:
    L = pickle.load(Lf)

with open(Trfile, 'rb') as Trf:
    Tr = pickle.load(Trf)

Tr.GMM_sample(1000000, L, simple=True)

figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)
Tr.gmm_samples = Tr.gmm_samples[:,0:9]

print "figure model len: ", len(figure_model.axes)

Tr.plot_kinematics('gmm', L.samples_means, L.samples_stds, color='m', fig=figure_model)
red_patch = mpatches.Patch(color='red')
magenta_patch = mpatches.Patch(color='magenta')
figure_model.legend(handles=[magenta_patch, red_patch], labels=["GMM Prior", "GMM Likelihood"], fontsize=16)

figure_model.savefig(save_path+'LvsPGMM.png', dpi=600)
figure_model.savefig(save_path+'LvsPGMM.pdf', dpi=600)

