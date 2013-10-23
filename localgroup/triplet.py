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
        Nsamples = 1
        
#         self.MW = localgroup.Halo('MW')
#         self.MW.sample_from(obs.data['MW'],Nsamples)
        
        self.M31 = localgroup.Halo('M31')
        self.M31.sample_from(obs.data['M31'],Nsamples)
        
        self.M33 = localgroup.Halo('M33')
        self.M33.sample_from(obs.data['M33'],Nsamples)
        
        return
        
# ----------------------------------------------------------------------------

    def get_parameters(self):
        
        return
        
# ----------------------------------------------------------------------------

    def get_properties(self):
        
        return
        
# ======================================================================

if __name__ == '__main__':

    t = Triplet()
    t.observe_halos()

    print "Checking M31:"
    # Check distances match:
    print t.M31.D,numpy.sqrt(t.M31.x*t.M31.x+t.M31.y*t.M31.y+t.M31.z*t.M31.z)
    # Check speeds match:
    u = numpy.sqrt(t.M31.v_west*t.M31.v_west + t.M31.v_north*t.M31.v_north + t.M31.vr*t.M31.vr)
    w = numpy.sqrt(t.M31.vx*t.M31.vx + t.M31.vy*t.M31.vy + t.M31.vz*t.M31.vz)
    print u,w

    print "Checking M33:"
    # Check distances match:
    print t.M33.D,numpy.sqrt(t.M33.x*t.M33.x+t.M33.y*t.M33.y+t.M33.z*t.M33.z)
    # Check speeds match:
    u = numpy.sqrt(t.M33.v_west*t.M33.v_west + t.M33.v_north*t.M33.v_north + t.M33.vr*t.M33.vr)
    w = numpy.sqrt(t.M33.vx*t.M33.vx + t.M33.vy*t.M33.vy + t.M33.vz*t.M33.vz)
    print u,w
