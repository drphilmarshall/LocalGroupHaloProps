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


figure_obs = L.plot_samples(10, color='b', overlay=False)

figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)

# In[6]:

L.model_gof(L.T.Nsamples, color="r", fig=figure_obs)
red_patch = mpatches.Patch(color='red')
blue_patch = mpatches.Patch(color='blue')
figure_obs.legend(handles=[blue_patch, red_patch], labels=["Observation Likelihood", "GMM Likelihood"])

figure_obs.savefig(save_path+"LGMM.pdf", dpi=800)
figure_obs.savefig(save_path+"LGMM.png", dpi=800)

Tr.GMM_sample(1000000, L, simple=True)


gmm_MW = np.copy(Tr.gmm_samples[:,10])
gmm_M31 = np.copy(Tr.gmm_samples[:,9])
gmm_M33 = np.copy(Tr.gmm_samples[:,11])
gmm_LMC = np.copy(Tr.gmm_samples[:,12])

gmm_MW_C = np.copy(Tr.gmm_samples[:,14])
gmm_M31_C = np.copy(Tr.gmm_samples[:,13])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))
Tr.gmm_samples = Tr.gmm_samples[:,0:9]


Tr.unprocess(L.samples_means, L.samples_stds, 'gmm')
data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_MW_C, gmm_M31_C)))
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{\\rm MW}$", "$Mvir_{\\rm M31}$", "$Mvir_{\\rm M33}$", "$Mvir_{\\rm LMC}$", "$Cvir_{\\rm MW}$", "$Cvir_{\\rm M31}$"]
pl = triangle.corner(data2, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                          plot_datapoints=False, bins=20, color='m')
Tr.preprocess(L.samples_means, L.samples_stds, mode='gmm')


Tr.unprocess(L.samples_means, L.samples_stds, mode='sim')
data = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.MW.Mvir), np.log10(Tr.M31.Mvir), np.log10(Tr.M33.Mvir), np.log10(Tr.LMC.Mvir), Tr.MW.Cvir, Tr.M31.Cvir)))
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{\\rm MW}$", "$Mvir_{\\rm M31}$", "$Mvir_{\\rm M33}$", "$Mvir_{\\rm LMC}$", "$Cvir_{\\rm MW}$", "$Cvir_{\\rm M31}$"]
sim_plot = triangle.corner(data, labels=labs, quantiles=[0.16,0.5,0.84], fig=pl, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 12},                          plot_datapoints=False, bins=20, color='c', label_kwargs={"fontsize": 16})
magenta_patch = mpatches.Patch(color='m')
cyan_patch = mpatches.Patch(color='c')
sim_plot.legend(handles=[cyan_patch, magenta_patch], labels=["CONSUELO Prior", "GMM-fit CONSUELO Prior"], fontsize=48)
Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')

sim_plot.savefig(save_path+'Q_GMMP_GOF.png', dpi=600)
sim_plot.savefig(save_path+'Q_GMMP_GOF.pdf', dpi=600)


