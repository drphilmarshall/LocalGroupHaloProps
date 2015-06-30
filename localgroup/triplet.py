# ======================================================================
from __future__ import division
import timingargument
import localgroup
import triangle
import numpy as np
import sys
import pickle
from sklearn import mixture
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Make sure Yao-Yuan Mao's "helpers" module is on your PYTHONPATH:
#   git clone git@bitbucket.org:yymao/helpers.git
from helpers.SimulationAnalysis import readHlist


Dvir_to_D200 = 0.047107689
D200_to_Dvir = 2.12279568
# ======================================================================

class Triplet(object):
    """
    NAME
        Triplet

    PURPOSE
        Define a local group analog, containing three halos, and
        perform various calculations.

    COMMENTS
        This class is used in two places: derivation of the local group
        analog likelihood function, and then inference of local group
        halo properties using that likelihood.

    INITIALISATION


    METHODS

        read_halos
        observe_halos
        get_parameters
        get_properties

    BUGS

    AUTHORS
      This file is part of the LocalGroupHaloProps project,
      distributed under the GPL v2,
      by Phil Marshall and Michael Busha (KIPAC).
      Please cite: Marshall et al in preparation.

    HISTORY
      2013-10-01  started Marshall & Busha (KIPAC)
    """
# ======================================================================

    def __init__(self, isPair=False):

        self.Nsamples = None
        self.isPair = isPair
        return

# ----------------------------------------------------------------------------

    def read_halos(self,directory):

        # Could be a lot more elegant.

        self.MW = localgroup.Halo()
        MWfile = directory+'/box01_MW.fits'
        self.MW.read(MWfile)

        self.M31 = localgroup.Halo()
        M31file = directory+'/box01_M31.fits'
        self.M31.read(MWfile)
        if not self.isPair:
            self.M33 = localgroup.Halo()
            M33file = directory+'/box01_M33.fits'
            self.M33.read(MWfile)

        return

# ----------------------------------------------------------------------------

    def observe_halos(self,Nsamples=10000):

        obs = localgroup.Observations()

        self.MW = localgroup.Halo('MW')
        self.MW.sample_from(obs.data['MW'],Nsamples)

        self.M31 = localgroup.Halo('M31')
        self.M31.sample_from(obs.data['M31'],Nsamples)
        if not self.isPair:
            self.M33 = localgroup.Halo('M33')
            self.M33.sample_from(obs.data['M33'],Nsamples)

            self.LMC = localgroup.Halo('LMC')
            self.LMC.sample_from(obs.data['LMC'], Nsamples)

        self.Nsamples = Nsamples

        return

# ----------------------------------------------------------------------------

# Transform from heliocentric spherical coordinates to M31-centric
# cartesian coordinates.

    def transform_to_M31(self, sim=False):

        # Covert M31 and M33 from heliocentric spherical to
        # galactocentric cartesian.
        if not sim:
            self.M31.x, self.M31.y, self.M31.z, self.M31.vx, self.M31.vy, self.M31.vz = localgroup.heliocentric_equatorial_spherical_to_galactocentric_cartesian(self.M31.RA,self.M31.DEC,self.M31.D,self.M31.v_west,self.M31.v_north,self.M31.v_r, self.M31.deltavrot_west, self.M31.deltavrot_north, R0=self.MW.x, VX=self.MW.vx, V0=self.MW.vy, VZ=self.MW.vz, M31=True)
            self.M31.frame = 'MW'
            if not self.isPair:
                self.M33.x, self.M33.y, self.M33.z, self.M33.vx, self.M33.vy, self.M33.vz = localgroup.heliocentric_equatorial_spherical_to_galactocentric_cartesian(self.M33.RA,self.M33.DEC,self.M33.D,self.M33.mu_west,self.M33.mu_north,self.M33.v_r, self.M33.deltavrot_west, self.M33.deltavrot_north, R0=self.MW.x, VX=self.MW.vx, V0=self.MW.vy, VZ=self.MW.vz, M31=False)
                self.M33.frame = 'MW'
                                
                self.LMC.x, self.LMC.y, self.LMC.z, self.LMC.vx, self.LMC.vy, self.LMC.vz = localgroup.heliocentric_equatorial_spherical_to_galactocentric_cartesian(self.LMC.RA,self.LMC.DEC,self.LMC.D,self.LMC.mu_west,self.LMC.mu_north,self.LMC.v_r, self.LMC.deltavrot_west, self.LMC.deltavrot_north, R0=self.MW.x, VX=self.MW.vx, V0=self.MW.vy, VZ=self.MW.vz, M31=False)
                self.LMC.frame = 'MW'
            # First we translate the MW positions from heliocentric
            # cartesian to galactocentric cartesian.
            #self.LMC.translate_to(self.MW)
            self.MW.translate_to(self.MW) #NOTE: This must be after heliocentric_equatorial_spherical_to_galactocentric_cartesian calls
            #self.LMC.translate_to(self.MW)
        # Now we can finally translate to M31 frame
        self.MW.translate_to(self.M31)
        if not self.isPair:
            self.M33.translate_to(self.M31)
            self.LMC.translate_to(self.M31)
        self.M31.translate_to(self.M31) #NOTE: This must be last
        if sim: self.sim_samples = np.transpose(np.array(self.get_kinematics()))
        return

# ============================================================================

    def read_sim_points(self, path, n_points, halo_props, h=0.7, a=1.0, npy=False, M31_larger=True):
        if npy:
            sim_data = np.load(path)
        else:
            sim_data = readHlist(path)
        if sim_data.shape[0] < n_points:
            raise ValueError('n_points too large.')
        np.random.shuffle(sim_data)
        self.sim_data = sim_data[:n_points]
        self.observe_halos(Nsamples=n_points)
        
        larger = 'M31_'
        smaller = 'MW_'
        larger_sub = 'M33_' #sub of larger host
        smaller_sub = 'LMC_' # sub of smaller host
        if not M31_larger:
            # Assign MW label to the larger halo in the data file
            larger = 'MW_'
            smaller = 'M31_'
            larger_sub = 'LMC_'
            smaller_sub = 'M33_'
        
        self.MW.translate_to(self.MW)
        self.MW.Mvir = h*self.sim_data[smaller+'mvir']
        self.MW.Rvir = h*self.sim_data[smaller+'rvir']
        self.MW.Rs = h*self.sim_data[smaller+'rs']
        self.MW.Cvir = self.MW.Rvir/self.MW.Rs

        self.M31.x = a*h*(self.sim_data[larger+'x'] - self.sim_data[smaller+'x'])
        self.M31.y = a*h*(self.sim_data[larger+'y'] - self.sim_data[smaller+'y'])
        self.M31.z = a*h*(self.sim_data[larger+'z'] - self.sim_data[smaller+'z'])
        self.M31.vx = self.sim_data[larger+'vx'] - self.sim_data[smaller+'vx']
        self.M31.vy = self.sim_data[larger+'vy'] - self.sim_data[smaller+'vy']
        self.M31.vz = self.sim_data[larger+'vz'] - self.sim_data[smaller+'vz']
        self.M31.frame = 'MW'
        self.M31.Mvir = h*self.sim_data[larger+'mvir']
        self.M31.Rvir = h*self.sim_data[larger+'rvir']
        self.M31.Rs = h*self.sim_data[larger+'rs']
        self.M31.Cvir = self.M31.Rvir/self.M31.Rs
        if not self.isPair:
            self.M33.x = a*h*(self.sim_data[larger_sub+'x'] - self.sim_data[smaller+'x'])
            self.M33.y = a*h*(self.sim_data[larger_sub+'y'] - self.sim_data[smaller+'y'])
            self.M33.z = a*h*(self.sim_data[larger_sub+'z'] - self.sim_data[smaller+'z'])
            self.M33.vx = self.sim_data[larger_sub+'vx'] - self.sim_data[smaller+'vx']
            self.M33.vy = self.sim_data[larger_sub+'vy'] - self.sim_data[smaller+'vy']
            self.M33.vz = self.sim_data[larger_sub+'vz'] - self.sim_data[smaller+'vz']
            self.M33.frame = 'MW'
            self.M33.Mvir = h*self.sim_data[larger_sub+'mvir']
            self.M33.Rvir = h*self.sim_data[larger_sub+'rvir']
            self.M33.Rs = h*self.sim_data[larger_sub+'rs']
            self.M33.Cvir = self.M33.Rvir/self.M33.Rs

            
            self.LMC.x = a*h*(self.sim_data[smaller_sub+'x'] - self.sim_data[smaller+'x'])
            self.LMC.y = a*h*(self.sim_data[smaller_sub+'y'] - self.sim_data[smaller+'y'])
            self.LMC.z = a*h*(self.sim_data[smaller_sub+'z'] - self.sim_data[smaller+'z'])
            self.LMC.vx = self.sim_data[smaller_sub+'vx'] - self.sim_data[smaller+'vx']
            self.LMC.vy = self.sim_data[smaller_sub+'vy'] - self.sim_data[smaller+'vy']
            self.LMC.vz = self.sim_data[smaller_sub+'vz'] - self.sim_data[smaller+'vz']
            self.LMC.frame = 'MW'
            self.LMC.Mvir = h*self.sim_data[smaller_sub+'mvir']
            self.LMC.Rvir = h*self.sim_data[smaller_sub+'rvir']
            self.LMC.Rs = h*self.sim_data[smaller_sub+'rs']
            self.LMC.Cvir = self.M33.Rvir/self.M33.Rs

            
        self.LG_Mvir = self.M31.Mvir + self.MW.Mvir
        return
# ============================================================================

    def mass_filter(self, mode='sim'):
        cond = self.MW.Mvir < self.M31.Mvir
        if mode == 'sim':
            self.sim_samples = self.sim_samples[cond]
            self.MW.Mvir = self.MW.Mvir[cond]
            self.M31.Mvir = self.M31.Mvir[cond]
            if not self.isPair: self.M33.Mvir = self.M33.Mvir[cond]
            self.LG_Mvir = self.LG_Mvir[cond]
            print 'sim_samples new shape: ', self.sim_samples.shape
        else:
            raise ValueError("Can only filter in mode=sim")
        return
# ============================================================================
    def dist_filter(self, condition):
        if not self.isPair:
            fig = plt.subplot(2,2,1)
            plt.hist(self.sim_samples[:,0])
            plt.title("MW Consuelo D distribution")
            plt.subplot(2,2,2)
            plt.hist(self.sim_samples[:,3])
            plt.title("M33 Consuelo D distribution")
            print "sim_sample length before: ", self.sim_samples.shape
            self.MW.Mvir = self.MW.Mvir[condition]
            self.MW.Rvir = self.MW.Rvir[condition]
            self.MW.Rs = self.MW.Rs[condition]
            self.MW.Cvir = self.MW.Cvir[condition]
            self.M31.Mvir = self.M31.Mvir[condition]
            self.M31.Rvir = self.M31.Rvir[condition]
            self.M31.Rs = self.M31.Rs[condition]
            self.M31.Cvir = self.M31.Cvir[condition]
            self.M33.Mvir = self.M33.Mvir[condition]
            self.M33.Rvir = self.M33.Rvir[condition]
            self.M33.Rs = self.M33.Rs[condition]
            self.M33.Cvir = self.M33.Cvir[condition]
            self.sim_samples = self.sim_samples[condition]
            print "sim_sample length after: ", self.sim_samples.shape
            plt.subplot(2,2,3)
            plt.hist(self.sim_samples[:,0])
            plt.title("MW Consuelo D distribution")
            plt.subplot(2,2,4)
            plt.hist(self.sim_samples[:,3])
            plt.title("M33 Consuelo D distribution")
        else:
            fig = plt.subplot(2,1,1)
            plt.hist(self.sim_samples[:,0])
            plt.title("MW Consuelo D distribution")
            print "sim_sample length before: ", self.sim_samples.shape
            self.MW.Mvir = self.MW.Mvir[condition]
            self.MW.Rvir = self.MW.Rvir[condition]
            self.MW.Rs = self.MW.Rs[condition]
            self.MW.Cvir = self.MW.Cvir[condition]
            self.M31.Mvir = self.M31.Mvir[condition]
            self.M31.Rvir = self.M31.Rvir[condition]
            self.M31.Rs = self.M31.Rs[condition]
            self.M31.Cvir = self.M31.Cvir[condition]
            self.sim_samples = self.sim_samples[condition]
            print "sim_sample length after: ", self.sim_samples.shape
            plt.subplot(2,1,2)
            plt.hist(self.sim_samples[:,0])
            plt.title("MW Consuelo D distribution")
        return fig

# ============================================================================

    def GMM(self, ngauss, data):
        self.gmm = mixture.GMM(ngauss, covariance_type='full')
        self.gmm_data = data
       # self.gmm_data_means = np.array([np.mean(self.gmm_data[:,i]) for i in range(self.gmm_data.shape[1])])
        #self.gmm_data_stds = np.array([np.std(self.gmm_data[:,i]) for i in range(self.gmm_data.shape[1])])
        #self.preprocess(self.gmm_data_means, self.gmm_data_stds, 'gmm_data')
        self.gmm.fit(data)
        return
# ============================================================================

    def GMM_sample(self, N):
        self.gmm_samples = self.gmm.sample(N)
        return
# ============================================================================

    def GMM_prior(self, ncomp, sample_size):
        dat = np.transpose(np.vstack((np.transpose(self.sim_samples), self.M31.Mvir, self.MW.Mvir, self.M33.Mvir)))
        self.GMM(ncomp, dat)
        self.GMM_sample(sample_size)
        self.unprocess(self.gmm_data_means, self.gmm_data_stds, 'gmm_data')
        self.unprocess(self.gmm_data_means, self.gmm_data_stds, 'gmm')
        self.MW.Mvir = np.copy(self.gmm_samples[:,7])
        self.M31.Mvir = np.copy(self.gmm_samples[:,6])
        self.LG_Mvir = np.copy(self.MW.Mvir) + np.copy(self.M31.Mvir)
        if not self.isPair:
            self.M33.Mvir = np.copy(self.gmm_samples[:,8])
        filter_nan = np.logical_not(np.isnan(self.MW.Mvir))&\
                     np.logical_not(np.isnan(self.M31.Mvir))&\
                     np.logical_not(np.isnan(self.LG_Mvir))
        if not self.isPair:
            filter_nan = filter_nan&np.logical_not(np.isnan(self.M33.Mvir))
        self.MW.Mvir = self.MW.Mvir[filter_nan]
        self.M31.Mvir = self.M31.Mvir[filter_nan]
        self.LG_Mvir = self.LG_Mvir[filter_nan]
        if not self.isPair:
            self.M33.Mvir = self.M33.Mvir[filter_nan]
        self.gmm_samples = self.gmm_samples[filter_nan]
        self.gmm_samples = self.gmm_samples[:,0:6]
        return
# ============================================================================

    def compute_timing_mass(self):
        M, a, x, e = timingargument.mass(self.MW.D, self.MW.v_r, vt=None, approach='radial', t0scatter=False)
        self.timing_masses = M
        return
# ============================================================================

    def compute_model_weights(self, L, mode, normalize=True):
        if mode == 'sim':
            weights = np.exp(L.evaluate(self.sim_samples)[0])
        elif mode == 'gmm':
            weights = np.exp(L.evaluate(self.gmm_samples)[0])
        #minw = weights.min()
        #weights = weights - minw
        total_weight = weights.sum()
        self.weights = 1.0/total_weight*weights
        return

# ============================================================================
    def calculate_N95(self):
        weights_copy = np.copy(self.weights)
        weights_copy.sort()
        weights_copy = weights_copy[::-1]
        sum = 0
        count = 0
        while sum < 0.95:
            sum = sum + weights_copy[count]
            count = count + 1
        return count
# ============================================================================
    def preprocess(self, means, stds, mode):
        if mode == 'sim':
            self.sim_samples = (self.sim_samples - means)/stds
        elif mode == 'gmm':
            self.gmm_samples = (self.gmm_samples - means)/stds
        elif mode == 'gmm_data':
            self.gmm_data = (self.gmm_data - means)/stds
        return

# ============================================================================

    def unprocess(self, means, stds, mode):
        if mode == 'sim':
            self.sim_samples = self.sim_samples*stds + means
        elif mode == 'gmm':
            self.gmm_samples = self.gmm_samples*stds + means
        elif mode == 'gmm_data':
            self.gmm_data = self.gmm_data*stds + means
        return

# ============================================================================

    def plot_kinematics(self, mode, means, stds, color, fig=None):
        self.unprocess(means, stds, mode)
        if mode == 'sim':
            data = self.sim_samples
        elif mode == 'gmm':
            data = self.gmm_samples[:,0:6]
        if self.isPair:
            # labs = ["MW_D", "MW_vr", "MW_vt"]
            labels = ["$D^{\\rm MW} Mpc$", "$v_{\\rm rad}^{\\rm MW} km/s$", "$v_{\\rm tan}^{\\rm MW} km/s$"]
        else:
            # labs = ["MW_D", "MW_vr", "MW_vt", "M33_D", "M33_vr", "M33_vt"]
            labels = ["$D^{\\rm MW} Mpc$", "$v_{\\rm rad}^{\\rm MW} km/s$", "$v_{\\rm tan}^{\\rm MW} km/s$", "$D^{\\rm M33} Mpc$", "$v_{\\rm rad}^{\\rm M33} km/s$", "$v_{\\rm tan}^{\\rm M33} km/s$"]

        if self.isPair:
            figure = triangle.corner(data, labels=labels, quantiles=[0.16,0.5,0.84], fig=fig, show_titles=True, title_args={"fontsize": 12}, label_args={"fontsize": 16}, color=color)
        else:
            figure = triangle.corner(data, labels=labels, quantiles=[0.16,0.5,0.84], fig=fig, show_titles=True, title_args={"fontsize": 12}, label_args={"fontsize": 16}, color=color)
        self.preprocess(means, stds, mode)
        return figure

# ============================================================================

    def save(self, save_path):
        save_file = open(save_path, "wb")
        pickle.dump(self, save_file)
        return

# ----------------------------------------------------------------------------

    def get_kinematics(self):
        if self.isPair:
            return self.MW.D, self.MW.v_r, self.MW.v_t
        else:
            return self.MW.D, self.MW.v_r, self.MW.v_t, self.M33.D, self.M33.v_r, self.M33.v_t, self.LMC.D, self.LMC.v_r, self.LMC.v_t

# ----------------------------------------------------------------------------

# Theory and expressions from Appendix in http://arxiv.org/pdf/0709.1159v1.pdf

    def calculate_cvir_to_c200(self, cvir):
        f = lambda x, y=0: x**3*(np.log(1.0+1.0/x) - 1.0/(1+x))-y
        inv_args = D200_to_Dvir*np.array(map(f, 1.0/cvir))
        c200 = 1.0/np.array([fsolve(f, [10], y)[0] for y in inv_args])
        return c200

    def calculate_Mvir_to_M200(self, Mvir, cvir, c200):
        M200 = Mvir*D200_to_Dvir*(c200/cvir)**3
        return np.array(M200)

# ============================================================================












if __name__ == '__main__':

    t = Triplet()
    t.observe_halos(Nsamples=200000)
    t.transform_to_M31()
    D_MW, vr_MW, vt_MW, D_M33, vr_M33, vt_M33, D_LMC, vr_LMC, vt_LMC = t.get_kinematics()

    print "Kinematics of first and 100th MW object: "
    i = 0
    print D_MW[i], vr_MW[i], vt_MW[i]
    i = 99
    print D_MW[i], vr_MW[i], vt_MW[i]

    print "Median MW kinematic parameters:"
    print np.median(D_MW), np.median(vr_MW), np.median(vt_MW)

    print "Median M33 kinematic parameters:"
    print np.median(D_M33), np.median(vr_M33), np.median(vt_M33)








    """

    print "Checking MW:"
    # Check speed:
    w = np.sqrt(t.MW.vx[0]*t.MW.vx[0] + t.MW.vy[0]*t.MW.vy[0] + t.MW.vz[0]*t.MW.vz[0])
    print "MW speed: ",w

    print "Checking M31:"
    # Check distances match:
    print "M31 distances: ",t.M31.D[0],np.sqrt(t.M31.x[0]*t.M31.x[0]+t.M31.y[0]*t.M31.y[0]+t.M31.z[0]*t.M31.z[0])
    # Check speeds match:
    u = np.sqrt(t.M31.v_west[0]*t.M31.v_west[0] + t.M31.v_north[0]*t.M31.v_north[0] + t.M31.v_r[0]*t.M31.v_r[0])
    w = np.sqrt(t.M31.vx[0]*t.M31.vx[0] + t.M31.vy[0]*t.M31.vy[0] + t.M31.vz[0]*t.M31.vz[0])
    print "M31 speeds: ",u,w

    print "Checking M33:"
    # Check distances match:
    print "M33 distances: ",t.M33.D[0],np.sqrt(t.M33.x[0]*t.M33.x[0]+t.M33.y[0]*t.M33.y[0]+t.M33.z[0]*t.M33.z[0])
    # Check speeds match:
    u = np.sqrt(t.M33.v_west[0]*t.M33.v_west[0] + t.M33.v_north[0]*t.M33.v_north[0] + t.M33.v_r[0]*t.M33.v_r[0])
    w = np.sqrt(t.M33.vx[0]*t.M33.vx[0] + t.M33.vy[0]*t.M33.vy[0] + t.M33.vz[0]*t.M33.vz[0])
    print "M33 speeds: ",u,w

    print " "
    """


    print "Before taking out solar motion, M31 velocity: ",np.mean(t.M31.vx),'+/-',np.std(t.M31.vx),', ', \
                           np.mean(t.M31.vy),'+/-',np.std(t.M31.vy),', ', \
                           np.mean(t.M31.vz),'+/-',np.std(t.M31.vz)

    print " "

    # Now transform to galactocentric coordinates:

    print "Calculating Galactocentric quantities..."
    #print np.mean(t.MW.vx), np.mean(t.MW.vy), np.mean(t.MW.vz)
    t.M31.translate_to(t.MW)
    t.M33.translate_to(t.MW)
    t.LMC.translate_to(t.MW)
    t.MW.translate_to(t.MW)

    # What is the space motion of M31? (eq 3 of vdM12)

    print "M31 position: ",np.mean(t.M31.x),np.mean(t.M31.y),np.mean(t.M31.z)
    print "  cf vdM++12: (-0.379, 0.613, -0.283)"

 #   print "M31 proper motion: ",np.mean(t.M31.v_west),'+/-',np.std(t.M31.v_west),', ', \
  #                              np.mean(t.M31.v_north),'+/-',np.std(t.M31.v_north)
   # print "  cf vdM++12: (-125+/-31, -74+/-28) km/s"
   # print "This is a little off because of our deltavrot hackery..."

    print " "

    print "M31 velocity: ",np.mean(t.M31.vx),'+/-',np.std(t.M31.vx),', ', \
                           np.mean(t.M31.vy),'+/-',np.std(t.M31.vy),', ', \
                           np.mean(t.M31.vz),'+/-',np.std(t.M31.vz)
    print "  cf vdM++12: (66+/-27, -76+/-19, 45+/-27) km/s"
    print "max vx = ",np.max(t.M31.vx)
    w = np.sqrt(t.M31.vx*t.M31.vx + t.M31.vy*t.M31.vy + t.M31.vz*t.M31.vz)
    print "M31 speed: ",np.mean(w),'+/-',np.std(w)
    print "  cf vdM++12: (110.6 +/- 7.8) km/s"


    print "M33 position: ",np.mean(t.M33.x),np.mean(t.M33.y),np.mean(t.M33.z)
    print "M33 velocity: ",np.mean(t.M33.vx),'+/-',np.std(t.M33.vx),', ', \
                           np.mean(t.M33.vy),'+/-',np.std(t.M33.vy),', ', \
                           np.mean(t.M33.vz),'+/-',np.std(t.M33.vz)

    print "LMC position: ",np.mean(t.LMC.x),np.mean(t.LMC.y),np.mean(t.LMC.z)
    print "LMC velocity: ",np.mean(t.LMC.vx),'+/-',np.std(t.LMC.vx),', ', \
                           np.mean(t.LMC.vy),'+/-',np.std(t.LMC.vy),', ', \
                           np.mean(t.LMC.vz),'+/-',np.std(t.LMC.vz)
