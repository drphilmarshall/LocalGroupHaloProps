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

save_path = "/lustre/ki/pfs/mwillia1/LG_project/plots/"
Lfile = '/lustre/ki/pfs/mwillia1/LG_project/L_Q_presym_LMCobs.pickle'
Trfile = '/lustre/ki/pfs/mwillia1/LG_project/Tr_sym_N30_v80_presym.pickle'
prop_sample_file = '/lustre/ki/pfs/mwillia1/LG_project/proposal_samples.pickle'


with open(Lfile, 'rb') as Lf:
    L = pickle.load(Lf)

with open(Trfile, 'rb') as Trf:
    Tr = pickle.load(Trf)

Tr.read_proposal_input_data(prop_sample_file)
Tr.fit_proposal_input(3, L)
Tr.sample_proposal(2000000)


# GOF of gmm fit to proposal distribution

Tr.unprocess(L.samples_means, L.samples_stds, 'prop')
data2 = Tr.proposal_samples
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{\\rm MW}$", "$Mvir_{\\rm M31}$", "$Mvir_{\\rm M33}$", "$Mvir_{\\rm LMC}$", "$Cvir_{\\rm MW}$", "$Cvir_{\\rm M31}$"]
pl = triangle.corner(data2, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                          plot_datapoints=False, bins=20, color='m')
Tr.preprocess(L.samples_means, L.samples_stds, mode='prop')


Tr.unprocess(L.samples_means, L.samples_stds, mode='prop_input')
data = Tr.input_prop_data
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{\\rm MW}$", "$Mvir_{\\rm M31}$", "$Mvir_{\\rm M33}$", "$Mvir_{\\rm LMC}$", "$Cvir_{\\rm MW}$", "$Cvir_{\\rm M31}$"]
sim_plot = triangle.corner(data, labels=labs, quantiles=[0.16,0.5,0.84], fig=pl, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 12},                          plot_datapoints=False, bins=20, color='c', label_kwargs={"fontsize": 16})
magenta_patch = mpatches.Patch(color='m')
cyan_patch = mpatches.Patch(color='c')
sim_plot.legend(handles=[cyan_patch, magenta_patch], labels=["Proposal Data", "GMM-fit Proposal"], fontsize=48)
Tr.preprocess(L.samples_means, L.samples_stds, mode='prop_input')

sim_plot.savefig(save_path+'Prop_vs_GMMProp.pdf', dpi=600)



# Comparison of prior with proposal








"""
gmm_MW = np.copy(Tr.proposal_samples[:,10])
gmm_M31 = np.copy(Tr.proposal_samples[:,9])
gmm_M33 = np.copy(Tr.proposal_samples[:,11])
gmm_LMC = np.copy(Tr.proposal_samples[:,12])

gmm_MW_C = np.copy(Tr.proposal_samples[:,14])
gmm_M31_C = np.copy(Tr.proposal_samples[:,13])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))
#Tr.gmm_samples = Tr.gmm_samples[:,0:9]

Tr.compute_model_weights(L, 'prop', normalize=True, split=1, imp=True)
count, w = Tr.calculate_N95(level=0.95, filter_samples=False, imp=True)

print "N95 = ", count
print "smallest weight = ", w



hlist_path = '/afs/slac.stanford.edu/u/ki/mwillia1/4001hlist.npy'
hlist = np.load(hlist_path)
hlist=hlist[np.abs(np.log10(hlist['mvir'])-12)<1.0]
hlist = hlist[hlist['upid']==-1]
cvir = hlist['rvir']/hlist['rs']
mvir = np.log10(hlist['mvir'])
bins = np.arange(np.min(mvir), np.max(mvir), .01)
dat = np.vstack((mvir, cvir)).T
conc = [[el[1] for el in dat if ((el[0] > bins[i]) & (el[0] < bins[i+1]))] for i in range(len(bins)-1)]
conc_means = np.array([np.mean(np.array(c)) for c in conc])
conc_stds = np.array([np.std(np.array(c)) for c in conc])







# In[ ]:

labs = ["$Mvir_{MW}$", "$Cvir_{MW}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_MW_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')

faxes = np.reshape(figure.axes, (2,2))
ax=faxes[1,0]
ax.plot(bins[1:], conc_means, lw=2, color='k')
ax.fill_between(bins[1:], conc_means+conc_stds, conc_means-conc_stds, facecolor='k', alpha=0.2)


figure.savefig(save_path+'Q_GMMP_MW_MvsC.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_MW_MvsC.pdf', dpi=800)


labs = ["$Mvir_{M31}$", "$Cvir_{M31}$"]
all_mvir = np.transpose(np.vstack((gmm_M31, gmm_M31_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')

faxes = np.reshape(figure.axes, (2,2))
ax=faxes[1,0]
ax.plot(bins[1:], conc_means, lw=2, color='k')
ax.fill_between(bins[1:], conc_means+conc_stds, conc_means-conc_stds, facecolor='k', alpha=0.2)


figure.savefig(save_path+'Q_GMMP_M31_MvsC.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_M31_MvsC.pdf', dpi=800)
# In[31]:

labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Mvir_{MW+M31}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_LG)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'Q_GMMP_all_Mvir.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_all_Mvir.pdf', dpi=800)
"""
