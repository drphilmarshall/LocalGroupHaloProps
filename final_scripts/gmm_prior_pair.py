
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

exist = int(sys.argv[1])
if exist:
    save_path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/working_plots/pair_exist/02/'
#    path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy' 
    path = '/afs/slac.stanford.edu/u/ki/mwillia1/All_Boxes_dat_M31_larger.npy'
else:
    save_path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/working_plots/pair/02/'
    path = '/afs/slac.stanford.edu/u/ki/mwillia1/All_Boxes_dat_M31_larger.npy'
#    path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes/Pairs/All_Boxes_dat_M31_larger.npy'
# Inside the Likelihood object is a "triplet" object called T, which contains an array of sample local groups, each with kinematic parameters consistent with the observational data. Let's plot these kinematic parameters in a "triangle" figure, to show all their 1 and 2-D marginal distributions.

# In[2]:

L = localgroup.Likelihood(isPair=True)
L.generate(Nsamples=200000)


# In[3]:

L.set_PDF(mixture.GMM(n_components=10, covariance_type='full'))
L.approximate()


# In[4]:

figure_obs = L.plot_samples(10, color='b', overlay=False)


# The above plot shows a Gaussian Mixture model fitted Gaussians. The shaded regions show two standard deviations. The samples data has been preprocessed to zero the mean and scale by standard deviation. Since we are using the Gaussian Mixture Model to model the underlying PDF of the data, more components is always better. 

# # How to evaluate goodness of fit:

# Due to lack of a standard goodness of fit test for GMM's, the best we can do is graphically show that the model reproduces the data well. We proceed by drawing a set of points from the fitted model, where each point is a local group with (MW_D, MW_vr, MW_vt, M33_D, M33_vr, M33_vt). We then plot the 1D and 2D marginalizations of the drawn point set and show that the marginalizations match the marginalizations of the true data.

# In[5]:

figure_model = L.model_gof(L.T.Nsamples, color="r", fig=None)


# In[6]:

L.model_gof(L.T.Nsamples, color="r", fig=figure_obs)
red_patch = mpatches.Patch(color='red')
blue_patch = mpatches.Patch(color='blue')
figure_obs.legend(handles=[red_patch, blue_patch], labels=["Model Generated", "Observation Generated"])


# In[7]:

figure_obs


# The above plot shows that the points drawn from the model create a population that is very similar to the true data.

# In[8]:

figure_obs.savefig(save_path+"LGMM.pdf", dpi=800)
figure_obs.savefig(save_path+"LGMM.png", dpi=800)

# # Reading Simulation Points:

# Below we read the preconfigured files containing the Consuelo (soon to be Dark Sky) Local Group analogs into a Triplet object. We plot the marginalizations of the simulation data, which allows us to compare with the LG prior.

# In[9]:

#path_cut = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy'
#path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy'

#path = '/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/data_files/MW_M31_pairs.txt'
npoints = 200000
halo_props = ['MW_Mvir', 'M31_Mvir', 'M33_Mvir']


# In[10]:
if exist:
    Tr = localgroup.Triplet(isPair=False)
else:
    Tr = localgroup.Triplet(isPair=True)

Tr.read_sim_points(path, npoints, halo_props, h=0.7, a=1.0, npy=True)

#Tr_cut = localgroup.Triplet(isPair=False)
#Tr_cut.read_sim_points(path_cut, npoints, halo_props, h=0.7, a=1.0, npy=True)


# In[11]:

Tr.transform_to_M31(sim=True)

#Tr_cut.transform_to_M31(sim=True)


# In[12]:

#Tr.mass_filter('sim')


# In[13]:

if exist:
    Tr.dist_filter((Tr.sim_samples[:,0] < 1) & (Tr.sim_samples[:,3] < 0.4) & (Tr.sim_samples[:,6] < 1))
else:
    Tr.dist_filter((Tr.sim_samples[:,0] < 1))

if exist: Tr.sim_samples = Tr.sim_samples[:,0:3]

# In[14]:

Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')

#Tr_cut.preprocess(L.samples_means, L.samples_stds, mode='sim')


# In[ ]:

sim_plot = Tr.plot_kinematics('sim', L.samples_means, L.samples_stds, color='c', fig=None)


# In[28]:

#sim_plot = Tr.plot_kinematics('sim', L.samples_means, L.samples_stds, color='c', fig=None)
# In[ ]:
if exist:
    dat = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.M31.Mvir), np.log10(Tr.MW.Mvir), np.log10(Tr.M33.Mvir), np.log10(Tr.LMC.Mvir))))
else:
    dat = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.M31.Mvir), np.log10(Tr.MW.Mvir))))

Tr.GMM(5, dat)

if exist:
    Tr_save_path = "/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr_exist.pickle"
    Tr.save(Tr_save_path) 
else:
    Tr_save_path = "/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr_pair.pickle"
    Tr.save(Tr_save_path)
#dat = np.transpose(np.vstack((np.transpose(Tr_cut.sim_samples), np.log10(Tr_cut.M31.Mvir), np.log10(Tr_cut.MW.Mvir), np.log10(Tr_cut.M33.Mvir))))
#Tr_cut.GMM(20, dat)


# In[ ]:

Tr.GMM_sample(500000, L, reps=1, simple=True)

#Tr_cut.GMM_sample(20000000)



# In[ ]:

gmm_MW = np.copy(Tr.gmm_samples[:,4])
gmm_M31 = np.copy(Tr.gmm_samples[:,3])
gmm_LG = np.log10(np.power(10,gmm_MW) + np.power(10,gmm_M31)) 
if exist:
    gmm_LMC = np.copy(Tr.gmm_samples[:,6])
    gmm_M33 = np.copy(Tr.gmm_samples[:,5])

Tr.gmm_samples = Tr.gmm_samples[:,0:3]



Tr.compute_model_weights(L, 'gmm', normalize=True, split=1)

#Tr_cut.compute_model_weights(L, 'gmm')


# In[ ]:

count, smallest_weight = Tr.calculate_N95(filter_samples=False)
print "N95 ", count
#print "smallest weight ", smallest_weight
#cond = Tr.weights[:] > smallest_weight
#gmm_MW = gmm_MW[cond]
#gmm_M31 = gmm_M31[cond]
#gmm_LG = gmm_LG[cond]
#if exist:
#    gmm_LMC = gmm_LMC[cond]
#    gmm_M33 = gmm_M33[cond]
#
#Tr.gmm_samples = Tr.gmm_samples[cond]
#Tr.weights = Tr.weights[cond]


Tr.unprocess(L.samples_means, L.samples_stds, 'gmm')
#data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31, gmm_M33, gmm_MW_C, gmm_M31_C, gmm_M33_C)))
if exist:
    data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31, gmm_M33, gmm_LMC)))
    labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$"]
else:
    data2 = np.transpose(np.vstack((np.transpose(Tr.gmm_samples), gmm_MW, gmm_M31)))
    labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$"]

pl = triangle.corner(data2, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                          plot_datapoints=False, bins=20, color='b')
Tr.preprocess(L.samples_means, L.samples_stds, mode='gmm')
# In[ ]:

if exist:
    data = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.MW.Mvir), np.log10(Tr.M31.Mvir), np.log10(Tr.M33.Mvir), np.log10(Tr.LMC.Mvir))))
    labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$"]
else:
    data = np.transpose(np.vstack((np.transpose(Tr.sim_samples), np.log10(Tr.MW.Mvir), np.log10(Tr.M31.Mvir))))
    labs = ["$D^{\\rm M31} Mpc$", "$v_{\\rm rad}^{\\rm M31} km/s$", "$v_{\\rm tan}^{\\rm M31} km/s$","$Mvir_{MW}$", "$Mvir_{M31}$"]


Tr.unprocess(L.samples_means, L.samples_stds, mode='sim')
sim_plot = triangle.corner(data, labels=labs, quantiles=[0.16,0.5,0.84], fig=pl, weights=None,                         plot_contours=True, show_titles=True, title_args={"fontsize": 12},                          plot_datapoints=False, bins=20, color='r', label_kwargs={"fontsize": 16})
red_patch = mpatches.Patch(color='r')
cyan_patch = mpatches.Patch(color='b')
sim_plot.legend(handles=[red_patch, cyan_patch], labels=["CONSUELO Prior", "GMM-fit CONSUELO Prior"], fontsize=48)
Tr.preprocess(L.samples_means, L.samples_stds, mode='sim')


# In[29]:



# In[ ]:

#name = 'gmm_CONSUELO_prior.png'
sim_plot.savefig(save_path+'P_GMMP_GOF.pdf', dpi=600)
"""
if exist:
    labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{M33}$", "$Mvir_{LMC}$", "$Mvir_{LG}$"]
    all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_M33, gmm_LMC, gmm_LG)))
else:
    labs = ["$Mvir_{MW}$", "$Mvir_{M31}$", "$Mvir_{LG}$"]
    all_mvir = np.transpose(np.vstack((gmm_MW, gmm_M31, gmm_LG)))

figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='c')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'P_GMMP_all_Mvir.pdf', dpi=800)


# In[ ]:

labs = ["$Mvir_{MW}$", "$Cvir_{MW}$"]
all_mvir = np.transpose(np.vstack((gmm_MW, gmm_MW_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='c')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'Q_GMMP_MW_MvsC.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_MW_MvsC.pdf', dpi=800)


labs = ["$Mvir_{M31}$", "$Cvir_{M31}$"]
all_mvir = np.transpose(np.vstack((gmm_M31, gmm_M31_C)))
figure = triangle.corner(all_mvir, labels=labs, quantiles=[0.16,0.5,0.84], fig=None, weights=Tr.weights,                         plot_contours=True, show_titles=True, title_args={"fontsize": 16}, label_args={"fontsize": 16},                         plot_datapoints=False, bins=20, color='c')
#figure.suptitle("Weighted Mass Posterior PDF, GMM Prior", fontsize=16, horizontalalignment='left')
figure.savefig(save_path+'Q_GMMP_M31_MvsC.png', dpi=800)
figure.savefig(save_path+'Q_GMMP_M31_MvsC.pdf', dpi=800)
# In[31]:
"""
