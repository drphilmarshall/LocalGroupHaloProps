
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

save_path = "/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/testplot/"
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

#figure_obs.savefig(save_path+"LGMM.pdf", dpi=800)
#figure_obs.savefig(save_path+"LGMM.png", dpi=800)

# # Reading Simulation Points:

# Below we read the preconfigured files containing the Consuelo (soon to be Dark Sky) Local Group analogs into a Triplet object. We plot the marginalizations of the simulation data, which allows us to compare with the LG prior.

# In[9]:

#path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/All_Boxes_quad_dat_M31_larger.npy'
path = '/afs/slac.stanford.edu/u/ki/mwillia1/All_Boxes_quad_symmetric_dat.npy'
npoints = 400000
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
Tr.GMM(30, dat)

Tr_save_path = "/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/Tr_sym_N30_v80_presym.pickle"
Tr.save(Tr_save_path)
with open("/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/L_Q_presym.pickle", "wb") as Lfile:
    pickle.dump(L, Lfile)


#dat = np.transpose(np.vstack((np.transpose(Tr_cut.sim_samples), np.log10(Tr_cut.M31.Mvir), np.log10(Tr_cut.MW.Mvir), np.log10(Tr_cut.M33.Mvir))))
#Tr_cut.GMM(20, dat)


# In[ ]:
