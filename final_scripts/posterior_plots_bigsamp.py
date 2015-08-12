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
gmm_sample_file1 = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/GMMsample1.bin'
gmm_sample_file2 = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/GMMsample2.bin'

with open(Lfile, 'rb') as Lf:
    L = pickle.load(Lf)

with open(Trfile, 'rb') as Trf:
    Tr = pickle.load(Trf)

with open(gmm_sample_file1, 'rb') as f:
    samples1 = pickle.load(f)

with open(gmm_sample_file2, 'rb') as g:
    samples2 = pickle.load(g)
r1, c1 = samples1.shape
r2, c2 = samples2.shape
samples1.resize((r1+r2, c1))
samples1[r1:] = samples2
Tr.GMM_sample(1, L, reps=1, simple=True)
Tr.gmm_samples = samples1


gmm_MW = np.copy(Tr.gmm_samples[:,10])
gmm_M31 = np.copy(Tr.gmm_samples[:,9])
gmm_M33 = np.copy(Tr.gmm_samples[:,11])
gmm_LMC = np.copy(Tr.gmm_samples[:,12])

gmm_MW_C = np.copy(Tr.gmm_samples[:,14])
gmm_M31_C = np.copy(Tr.gmm_samples[:,13])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))
Tr.gmm_samples = Tr.gmm_samples[:,0:9]

Tr.compute_model_weights(L, 'gmm', normalize=True)



# begin N95 calculation and filtering

#count, smallest_weight = Tr.calculate_N95(filter_samples=False)
#print "N95 ", count
#print "smallest weight ", smallest_weight


#cond = Tr.weights[:] > smallest_weight
#gmm_MW = gmm_MW[cond]
#gmm_M31 = gmm_M31[cond]
#gmm_M33 = gmm_M33[cond]
#gmm_LMC = gmm_LMC[cond]
#gmm_LG = gmm_LG[cond]

#gmm_MW_C = gmm_MW_C[cond]
#gmm_M31_C = gmm_M31_C[cond]
#Tr.gmm_samples = Tr.gmm_samples[cond]
#Tr.weights = Tr.weights[cond]

# end N95 calculation and filtering

hlist_path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/4001/4001hlist.npy'
hlist = np.load(hlist_path)
hlist=hlist[np.abs(np.log10(hlist['mvir'])-12)<1.0]
cvir = hlist['rvir']/hlist['rs']
mvir = np.log10(hlist['mvir'])
bins = np.arange(np.min(mvir), np.max(mvir), .01)
dat = np.vstack((mvir, cvir)).T
conc = [[el[1] for el in dat if ((el[0] > bins[i]) & (el[0] < bins[i+1]))] for i in range(len(bins)-1)]
conc_means = np.array([np.mean(np.array(c)) for c in conc])
conc_stds = np.array([np.std(np.array(c)) for c in conc])


labs = ["$Mvir_{MW}$", "$Cvir_{MW}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_MW_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')

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

labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Mvir_{LG}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_LG)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'Q_GMMP_all_Mvir.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_all_Mvir.pdf', dpi=800)
