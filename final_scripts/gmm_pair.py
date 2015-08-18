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

save_path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/working_plots/pair/02/'

L = localgroup.Likelihood(isPair=True)
L.generate(Nsamples=200000)

L.set_PDF(mixture.GMM(n_components=10, covariance_type='full'))
L.approximate()

figure_obs = L.plot_samples(10, color='b', overlay=False)

figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)

L.model_gof(L.T.Nsamples, color="r", fig=figure_obs)
red_patch = mpatches.Patch(color='red')
blue_patch = mpatches.Patch(color='blue')
figure_obs.legend(handles=[red_patch, blue_patch], labels=["Model Generated", "Observation Generated"])

path = '/afs/slac.stanford.edu/u/ki/mwillia1/All_Boxes_dat_M31_larger.npy'
npoints = 200000
halo_props = ['MW_Mvir', 'M31_Mvir', 'M33_Mvir']

Tr = localgroup.Triplet(isPair=True)
Tr.read_sim_points(path, npoints, halo_props, h=0.7, a=1.0, npy=True)

Tr.transform_to_M31(sim=True)

Tr.dist_filter((Tr.sim_samples[:,0] < 10))

Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')

dat = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.M31.Mvir), np.log10(Tr.MW.Mvir))))
Tr.GMM(20, dat)

Tr.GMM_sample(10000000, L, reps=1, simple=False)

path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr_P_samp.pickle'
with open(save_path, 'wb') as f:
    pickle.dump(Tr.gmm_samples, f)

"""
gmm_MW = np.copy(Tr.gmm_samples[:,4])
gmm_M31 = np.copy(Tr.gmm_samples[:,3])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))

Tr.gmm_samples = Tr.gmm_samples[:,0:3]

Tr.compute_model_weights(L, 'gmm')

count, smallest_weight = Tr.calculate_N95(filter_samples=False)
print "N95 ", count

Tr.unprocess(L.samples_means, L.samples_stds, 'gmm')
data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31)))
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$"]
pl = triangle.corner(data2, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=None,\
                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},\
                         plot_datapoints=False, bins=50, color='b')
Tr.preprocess(L.samples_means, L.samples_stds, mode='gmm')

Tr.unprocess(L.samples_means, L.samples_stds, mode='sim')
data = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.MW.Mvir), np.log10(Tr.M31.Mvir))))
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$"]
sim_plot = triangle.corner(data, labels=labs, quantiles=[0.16,0.5,0.84], fig=pl, weights=None,\
                         plot_contours=True, show_titles=False, title_args={"fontsize": 12}, \
                         plot_datapoints=False, bins=50, color='r', label_kwargs={"fontsize": 16}, label_args={"fontsize": 16})
red_patch = mpatches.Patch(color='r')
cyan_patch = mpatches.Patch(color='b')
sim_plot.legend(handles=[red_patch, cyan_patch], labels=["CONSUELO Prior", "GMM-fit CONSUELO Prior"], fontsize=16)
sim_plot.suptitle("GMM Prior Overlayed with Consuelo Prior", fontsize=16)
Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')

sim_plot.savefig(save_path+'P_GMMP_GOF.pdf', dpi=600)




labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{LG}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_LG)))

figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='c')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'P_GMMP_all_Mvir.pdf', dpi=800)
"""
