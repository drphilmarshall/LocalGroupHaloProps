# ======================================================================

import localgroup

import numpy

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

    def __init__(self):

        self.Nsamples = 100000
    
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
        
        self.M33 = localgroup.Halo()
        M33file = directory+'/box01_M33.fits'
        self.M33.read(MWfile)
        
        return
        
# ----------------------------------------------------------------------------

    def observe_halos(self):
        
        obs = localgroup.Observations()
        
        self.MW = localgroup.Halo('MW')
        self.MW.sample_from(obs.data['MW'],self.Nsamples)
        
        self.M31 = localgroup.Halo('M31')
        self.M31.sample_from(obs.data['M31'],self.Nsamples)
        
        self.M33 = localgroup.Halo('M33')
        self.M33.sample_from(obs.data['M33'],self.Nsamples)
        
        return
        
# ----------------------------------------------------------------------------

# Transform from heliocentric spherical coordinates to M31-centric
# cartesian coordinates.  

    def transform_to_M31(self):

        # Covert M31 and M33 from heliocentric spherical to 
        # galactocentric cartesian.  
        self.M31.x, self.M31.y, self.M31.z, self.M31.vx, self.M31.vy, self.M31.vz = localgroup.heliocentric_equatorial_spherical_to_galactocentric_cartesian(self.M31.RA,self.M31.DEC,self.M31.D,self.M31.mu_west,self.M31.mu_north,self.M31.v_r, R0=self.MW.x, V0=self.MW.vy)
        self.M31.frame = 'MW'
        self.M33.x, self.M33.y, self.M33.z, self.M33.vx, self.M33.vy, self.M33.vz = localgroup.heliocentric_equatorial_spherical_to_galactocentric_cartesian(self.M33.RA,self.M33.DEC,self.M33.D,self.M33.mu_west,self.M33.mu_north,self.M33.v_r, R0=self.MW.x, V0=self.MW.vy)
        self.M33.frame = 'MW'

        # First we translate the MW positions from heliocentric
        # cartesian to galactocentric cartesian.  
        self.MW.translate_to(self.MW) #NOTE: This must be after heliocentric_equatorial_spherical_to_galactocentric_cartesian calls

        # Now we can finally translate to M31 frame
        self.MW.translate_to(self.M31)
        self.M33.translate_to(self.M31)
        self.M31.translate_to(self.M31) #NOTE: This must be last

        return

# ----------------------------------------------------------------------------

    def get_parameters(self):
        
        return self.MW.D, self.MW.v_r, self.MW.v_t, self.M33.D, self.M33.v_r, self.M33.v_t,
        
# ----------------------------------------------------------------------------

    def get_properties(self):
        
        return

# ----------------------------------------------------------------------------

    def set_Nsamples(self, Nsamples):

        self.Nsamples = Nsamples

        return

        
# ======================================================================

if __name__ == '__main__':

    t = Triplet()
    t.observe_halos()
    t.transform_to_M31()
    D_MW, vr_MW, vt_MW, D_M33, vr_M33, vt_M33 = t.get_parameters()

    print "Parameters of first and 100th MW object: "
    i = 0
    print D_MW[i], vr_MW[i], vt_MW[i]
    i = 99
    print D_MW[i], vr_MW[i], vt_MW[i]
 
    print "Median MW Parameters:"
    print numpy.median(D_MW), numpy.median(vr_MW), numpy.median(vt_MW)

    print "Median M33 Parameters:"
    print numpy.median(D_M33), numpy.median(vr_M33), numpy.median(vt_M33)


    """
    print "Checking MW:"
    # Check speed:
    w = numpy.sqrt(t.MW.vx[0]*t.MW.vx[0] + t.MW.vy[0]*t.MW.vy[0] + t.MW.vz[0]*t.MW.vz[0])
    print "MW speed: ",w

    print "Checking M31:"
    # Check distances match:
    print "M31 distances: ",t.M31.D[0],numpy.sqrt(t.M31.x[0]*t.M31.x[0]+t.M31.y[0]*t.M31.y[0]+t.M31.z[0]*t.M31.z[0])
    # Check speeds match:
    u = numpy.sqrt(t.M31.v_west[0]*t.M31.v_west[0] + t.M31.v_north[0]*t.M31.v_north[0] + t.M31.v_r[0]*t.M31.v_r[0])
    w = numpy.sqrt(t.M31.vx[0]*t.M31.vx[0] + t.M31.vy[0]*t.M31.vy[0] + t.M31.vz[0]*t.M31.vz[0])
    print "M31 speeds: ",u,w

    print "Checking M33:"
    # Check distances match:
    print "M33 distances: ",t.M33.D[0],numpy.sqrt(t.M33.x[0]*t.M33.x[0]+t.M33.y[0]*t.M33.y[0]+t.M33.z[0]*t.M33.z[0])
    # Check speeds match:
    u = numpy.sqrt(t.M33.v_west[0]*t.M33.v_west[0] + t.M33.v_north[0]*t.M33.v_north[0] + t.M33.v_r[0]*t.M33.v_r[0])
    w = numpy.sqrt(t.M33.vx[0]*t.M33.vx[0] + t.M33.vy[0]*t.M33.vy[0] + t.M33.vz[0]*t.M33.vz[0])
    print "M33 speeds: ",u,w

    print " "
    
    print "Before taking out solar motion, M31 velocity: ",numpy.mean(t.M31.vx),'+/-',numpy.std(t.M31.vx),', ', \
                           numpy.mean(t.M31.vy),'+/-',numpy.std(t.M31.vy),', ', \
                           numpy.mean(t.M31.vz),'+/-',numpy.std(t.M31.vz)
   
    print " "

    # Now transform to galactocentric coordinates:
    
    print "Calculating Galactocentric quantities..."
    t.M31.translate_to(t.MW)
    t.M33.translate_to(t.MW)
    t.MW.translate_to(t.MW)
    
    # What is the space motion of M31? (eq 3 of vdM12)
    
    print "M31 position: ",numpy.mean(t.M31.x),numpy.mean(t.M31.y),numpy.mean(t.M31.z)
    print "  cf vdM++12: (-0.379, 0.613, -0.283)"
    
    print "M31 proper motion: ",numpy.mean(t.M31.v_west),'+/-',numpy.std(t.M31.v_west),', ', \
                                numpy.mean(t.M31.v_north),'+/-',numpy.std(t.M31.v_north)
    print "  cf vdM++12: (-125+/-31, -74+/-28) km/s"
    print "This is a little off because of our deltavrot hackery..."
   
    print " "
    
    print "M31 velocity: ",numpy.mean(t.M31.vx),'+/-',numpy.std(t.M31.vx),', ', \
                           numpy.mean(t.M31.vy),'+/-',numpy.std(t.M31.vy),', ', \
                           numpy.mean(t.M31.vz),'+/-',numpy.std(t.M31.vz)
    print "  cf vdM++12: (66+/-27, -76+/-19, 45+/-27) km/s"
    
    w = numpy.sqrt(t.M31.vx*t.M31.vx + t.M31.vy*t.M31.vy + t.M31.vz*t.M31.vz)
    print "M31 speed: ",numpy.mean(w),'+/-',numpy.std(w)
    print "  cf vdM++12: (110.6 +/- 7.8) km/s"

    print "Looks like there's a bug in the Galactocentric velocity transformation..."
    """
