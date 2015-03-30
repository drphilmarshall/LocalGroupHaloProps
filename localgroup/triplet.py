# ======================================================================
import timingargument
import localgroup
import triangle
import numpy as np
import sys
import pickle

sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist
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

            # First we translate the MW positions from heliocentric
            # cartesian to galactocentric cartesian.
            self.MW.translate_to(self.MW) #NOTE: This must be after heliocentric_equatorial_spherical_to_galactocentric_cartesian calls

        # Now we can finally translate to M31 frame
        self.MW.translate_to(self.M31)
        if not self.isPair:
            self.M33.translate_to(self.M31)
        self.M31.translate_to(self.M31) #NOTE: This must be last
        if sim: self.sim_samples = np.transpose(np.array(self.get_kinematics()))
        return

# ============================================================================

    def read_sim_points(self, path, n_points, halo_props):

        sim_data = readHlist(path)
        if sim_data.shape[0] < n_points:
            raise ValueError('n_points too large.')
        np.random.shuffle(sim_data)
        self.sim_data = sim_data[:n_points]
        self.observe_halos(Nsamples=n_points)

        self.MW.translate_to(self.MW)
        self.MW.Mvir = self.sim_data['MW_Mvir']

        self.M31.x = self.sim_data['M31_x'] - self.sim_data['MW_x']
        self.M31.y = self.sim_data['M31_y'] - self.sim_data['MW_y']
        self.M31.z = self.sim_data['M31_z'] - self.sim_data['MW_z']
        self.M31.vx = self.sim_data['M31_vx'] - self.sim_data['MW_vx']
        self.M31.vy = self.sim_data['M31_vy'] - self.sim_data['MW_vy']
        self.M31.vz = self.sim_data['M31_vz'] - self.sim_data['MW_vz']
        self.M31.frame = 'MW'
        self.M31.Mvir = self.sim_data['M31_Mvir']
        if not self.isPair:
            self.M33.x = self.sim_data['M33_x'] - self.sim_data['MW_x']
            self.M33.y = self.sim_data['M33_y'] - self.sim_data['MW_y']
            self.M33.z = self.sim_data['M33_z'] - self.sim_data['MW_z']
            self.M33.vx = self.sim_data['M33_vx'] - self.sim_data['MW_vx']
            self.M33.vy = self.sim_data['M33_vy'] - self.sim_data['MW_vy']
            self.M33.vz = self.sim_data['M33_vz'] - self.sim_data['MW_vz']
            self.M33.frame = 'MW'
            self.M33.Mvir = self.sim_data['M33_Mvir']
        return
# ============================================================================

    def compute_timing_mass(self):
        M, a, x, e = timingargument.mass(self.MW.D, self.MW.v_r, vt=None, approach='radial', t0scatter=False)
        self.timing_masses = M
        return
# ============================================================================

    def compute_model_weights(self, L, normalize=True):
        weights = np.exp(L.evaluate(self.sim_samples)[0])
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
    def preprocess(self, means, stds):
        self.sim_samples = (self.sim_samples - means)/stds

        return
# ============================================================================

    def plot_kinematics(self):

        if self.isPair:
            # labs = ["MW_D", "MW_vr", "MW_vt"]
            labels = ["$\Delta D^{\\rm MW}$", "$\Delta v_{\\rm rad}^{\\rm MW}$", "$\Delta v_{\\rm tan}^{\\rm MW}$"]
        else:
            # labs = ["MW_D", "MW_vr", "MW_vt", "M33_D", "M33_vr", "M33_vt"]
            labels = ["$\Delta D^{\\rm MW}$", "$\Delta v_{\\rm rad}^{\\rm MW}$", "$\Delta v_{\\rm tan}^{\\rm MW}$", "$\Delta D^{\\rm M33}$", "$\Delta v_{\\rm rad}^{\\rm M33}$", "$\Delta v_{\\rm tan}^{\\rm M33}$"]

        if self.isPair:
            figure = triangle.corner(self.sim_samples, labels=labels, quantiles=[0.16,0.5,0.84], show_titles=True, title_args={"fontsize": 12})
        else:
            figure = triangle.corner(self.sim_samples, labels=labels, quantiles=[0.16,0.5,0.84], show_titles=True, title_args={"fontsize": 12})

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
            return self.MW.D, self.MW.v_r, self.MW.v_t, self.M33.D, self.M33.v_r, self.M33.v_t

# ----------------------------------------------------------------------------

if __name__ == '__main__':

    t = Triplet()
    t.observe_halos(Nsamples=200000)
    t.transform_to_M31()
    D_MW, vr_MW, vt_MW, D_M33, vr_M33, vt_M33 = t.get_kinematics()

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
