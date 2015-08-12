
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


# Inside the Likelihood object is a "triplet" object called T, which contains an array of sample local groups, each with kinematic parameters consistent with the observational data. Let's plot these kinematic parameters in a "triangle" figure, to show all their 1 and 2-D marginal distributions.

# In[2]:

L = localgroup.Likelihood(isPair=False)
L.generate(Nsamples=200000)


# In[3]:

L.set_PDF(mixture.GMM(n_components=10, covariance_type='full'))
L.approximate()


# In[4]:

figure_obs = L.plot_samples(10, color='m', overlay=False)


# The above plot shows a Gaussian Mixture model fitted Gaussians. The shaded regions show two standard deviations. The samples data has been preprocessed to zero the mean and scale by standard deviation. Since we are using the Gaussian Mixture Model to model the underlying PDF of the data, more components is always better. 

# # How to evaluate goodness of fit:

# Due to lack of a standard goodness of fit test for GMM's, the best we can do is graphically show that the model reproduces the data well. We proceed by drawing a set of points from the fitted model, where each point is a local group with (MW_D, MW_vr, MW_vt, M33_D, M33_vr, M33_vt). We then plot the 1D and 2D marginalizations of the drawn point set and show that the marginalizations match the marginalizations of the true data.

# In[5]:

figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)
print "made figure model: ", len(figure_model.axes)
print figure_model.axes
#figure_model.savefig("/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/figmodel.pdf", dpi=800)

# In[6]:

L.model_gof(L.T.Nsamples, color="r", fig=figure_obs)
red_patch = mpatches.Patch(color='red')
magenta_patch = mpatches.Patch(color='magenta')
figure_obs.legend(handles=[red_patch, magenta_patch], labels=["Model Generated", "Observation Generated"])


# In[7]:

figure_obs


# The above plot shows that the points drawn from the model create a population that is very similar to the true data.

# In[8]:

figure_obs.savefig("/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/LGMM.pdf", dpi=800)
figure_obs.savefig("/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/LGMM.png", dpi=800)

# # Reading Simulation Points:

# Below we read the preconfigured files containing the Consuelo (soon to be Dark Sky) Local Group analogs into a Triplet object. We plot the marginalizations of the simulation data, which allows us to compare with the LG prior.

# In[9]:

path_cut = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy'
path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy'
#path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/data_files/MW_M31_pairs.txt'
npoints = 200000
halo_props = ['MW_Mvir', 'M31_Mvir', 'M33_Mvir']


# In[10]:

Tr = localgroup.Triplet(isPair=False)
Tr.read_sim_points(path, npoints, halo_props, h=0.7, a=1.0, npy=True)

#Tr_cut = localgroup.Triplet(isPair=False)
#Tr_cut.read_sim_points(path_cut, npoints, halo_props, h=0.7, a=1.0, npy=True)


# In[11]:

Tr.transform_to_M31(sim=True)

#Tr_cut.transform_to_M31(sim=True)


# In[12]:

#Tr.mass_filter('sim')


# In[13]:

Tr.dist_filter((Tr.sim_samples[:,0] < 1) & (Tr.sim_samples[:,3] < 0.4) & (Tr.sim_samples[:,6] < 1))

#Tr_cut.dist_filter((Tr_cut.sim_samples[:,0] < 1) & (Tr_cut.sim_samples[:,3] < 0.4))


# In[14]:

Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')

#Tr_cut.preprocess(L.samples_means, L.samples_stds, mode='sim')


# In[ ]:

sim_plot = Tr.plot_kinematics('sim', L.samples_means, L.samples_stds, color='b', fig=None)


# In[28]:

#sim_plot = Tr.plot_kinematics('sim', L.samples_means, L.samples_stds, color='c', fig=None)
# In[ ]:
print Tr.M33.Mvir.shape
print Tr.LMC.Mvir.shape
dat = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.M31.Mvir), np.log10(Tr.MW.Mvir), np.log10(Tr.M33.Mvir), np.log10(Tr.LMC.Mvir), Tr.M31.Cvir, Tr.MW.Cvir)))
#dat = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.M31.Mvir), np.log10(Tr.MW.Mvir))))
Tr.GMM(20, dat)

#dat = np.transpose(np.vstack((np.transpose(Tr_cut.sim_samples), np.log10(Tr_cut.M31.Mvir), np.log10(Tr_cut.MW.Mvir), np.log10(Tr_cut.M33.Mvir))))
#Tr_cut.GMM(20, dat)


# In[ ]:

Tr.GMM_sample(1000000, L)


#Tr_cut.GMM_sample(20000000)



# In[ ]:

gmm_MW = np.copy(Tr.gmm_samples[:,10])
gmm_M31 = np.copy(Tr.gmm_samples[:,9])
gmm_M33 = np.copy(Tr.gmm_samples[:,11])
gmm_LMC = np.copy(Tr.gmm_samples[:,12])

gmm_MW_C = np.copy(Tr.gmm_samples[:,14])
gmm_M31_C = np.copy(Tr.gmm_samples[:,13])
#gmm_M33_C = np.copy(Tr.gmm_samples[:,11])

gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))

#cond = gmm_MW < gmm_M31
#Tr.gmm_samples = Tr.gmm_samples[cond]
#gmm_MW = gmm_MW[cond]
#gmm_M31 = gmm_M31[cond]
#gmm_M33 = gmm_M33[cond]
#gmm_LG = gmm_LG[cond]


figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)
#print "figmod: ", len(figmod.axes)





Tr.gmm_samples = Tr.gmm_samples[:,0:9]

print "figure model len: ", len(figure_model.axes)
#print figure_model.axes

Tr.plot_kinematics('gmm', L.samples_means, L.samples_stds, color='c', fig=figure_model)
red_patch = mpatches.Patch(color='red')
cyan_patch = mpatches.Patch(color='cyan')
figure_model.legend(handles=[cyan_patch, red_patch], labels=["GMM Prior", "GMM Likelihood"], fontsize=16)

figure_model.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/LvsPGMM.png', dpi=600)
figure_model.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/LvsPGMM.pdf', dpi=600)

Tr.GMM_sample(200000000, L, repeat=2)

gmm_MW = np.copy(Tr.gmm_samples[:,10])
gmm_M31 = np.copy(Tr.gmm_samples[:,9])
gmm_M33 = np.copy(Tr.gmm_samples[:,11])
gmm_LMC = np.copy(Tr.gmm_samples[:,12])

gmm_MW_C = np.copy(Tr.gmm_samples[:,14])
gmm_M31_C = np.copy(Tr.gmm_samples[:,13])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31))
Tr.gmm_samples = Tr.gmm_samples[:,0:9]

Tr.compute_model_weights(L, 'gmm', normalize=True, split=40)

#Tr_cut.compute_model_weights(L, 'gmm')


# In[ ]:

count, smallest_weight = Tr.calculate_N95(filter_samples=False)
print "N95 ", count
print "smallest weight ", smallest_weight
cond = Tr.weights[:] > smallest_weight
gmm_MW = gmm_MW[cond]
gmm_M31 = gmm_M31[cond]
gmm_M33 = gmm_M33[cond]
gmm_LMC = gmm_LMC[cond]
gmm_LG = gmm_LG[cond]

gmm_MW_C = gmm_MW_C[cond]
gmm_M31_C = gmm_M31_C[cond]
Tr.gmm_samples = Tr.gmm_samples[cond]
Tr.weights = Tr.weights[cond]
#print Tr_cut.calculate_N95()


# In[18]:


Tr.unprocess(L.samples_means, L.samples_stds, 'gmm')
#data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31, gmm_M33, gmm_MW_C, gmm_M31_C, gmm_M33_C)))
data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_MW_C, gmm_M31_C)))
#labs=["mwd", "mwvr", "mwvt", "m33d", "m33vr", "m33vt", "MWMvir", "M31Mvir", "M33Mvir"]
#labs=["mwd", "mwvr", "mwvt", "m33d", "m33vr", "m33vt", "LMCd", "LMCvr", "LMCvt", "MWMvir", "M31Mvir"]
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Cvir_{MW}$", "$Cvir_{M31}$"]
#labs=["mwd", "mwvr", "mwvt", "m33d", "m33vr", "m33vt", "MWMvir", "M31Mvir", "M33Mvir", "MW Cvir", "M31 Cvir", "M33 Cvir"]
pl = triangle.corner(data2, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                          plot_datapoints=False, bins=20, color='c')
Tr.preprocess(L.samples_means, L.samples_stds, mode='gmm')
# In[ ]:

Tr.unprocess(L.samples_means, L.samples_stds, mode='sim')
data = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.MW.Mvir), np.log10(Tr.M31.Mvir), np.log10(Tr.M33.Mvir), np.log10(Tr.LMC.Mvir), Tr.MW.Cvir, Tr.M31.Cvir)))
labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$","$D^{\\rm LMC} Mpc$", "$v_{\\rm rad}^{\\rm LMC} km/s$", "$v_{\\rm tan}^{\\rm LMC} km/s$", "$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Cvir_{MW}$", "$Cvir_{M31}$"]
sim_plot = triangle.corner(data, labels=labs, quantiles=[0.16,0.5,0.84], fig=pl, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 12},                          plot_datapoints=False, bins=20, color='b', label_kwargs={"fontsize": 16})
blue_patch = mpatches.Patch(color='b')
cyan_patch = mpatches.Patch(color='c')
sim_plot.legend(handles=[blue_patch, cyan_patch], labels=["CONSUELO Prior", "GMM-fit CONSUELO Prior"], fontsize=48)
Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')


# In[29]:



# In[ ]:

#name = 'gmm_CONSUELO_prior.png'
sim_plot.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_GOF.png', dpi=600)
sim_plot.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_GOF.pdf', dpi=600)



# In[ ]:

labs = ["$Mvir_{MW}$", "$Cvir_{MW}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_MW_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_MW_MvsC.png', dpi=800)
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_MW_MvsC.pdf', dpi=800)


labs = ["$Mvir_{M31}$", "$Cvir_{M31}$"]
all_mvir = np.transpose(np.vstack((gmm_M31, gmm_M31_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_M31_MvsC.png', dpi=800)
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_M31_MvsC.pdf', dpi=800)
# In[31]:

labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Mvir_{MW+M31}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_LG)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='g')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_all_Mvir.png', dpi=800)
figure.savefig('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Q_GMMP_all_Mvir.pdf', dpi=800)


