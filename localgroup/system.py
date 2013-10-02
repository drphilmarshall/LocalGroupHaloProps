# ======================================================================

import localgroup

import numpy

# ======================================================================

class System(object):
    """
    NAME
        System

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
        
        MW = localgroup.Halo()
        MWfile = directory+'/box01_MW.fits'
        MW.read(MWfile)
        
        M31 = localgroup.Halo()
        M31file = directory+'/box01_M31.fits'
        M31.read(MWfile)
        
        M33 = localgroup.Halo()
        M33file = directory+'/box01_M33.fits'
        M33.read(MWfile)
        
        return
        
# ----------------------------------------------------------------------------

    def observe_halos(self):
        
        obs = localgroup.Observations()
        Nsamples = 1000
        
        MW = localgroup.Halo()
        MW.sample_from(obs.data['MW'],Nsamples)
        
        M31 = localgroup.Halo()
        M31.sample_from(obs.data['M31'],Nsamples)
        
        M33 = localgroup.Halo()
        M33.sample_from(obs.data['M33'],Nsamples)
        
        return
        
# ----------------------------------------------------------------------------

    def get_parameters(self):
        
        return
        
# ----------------------------------------------------------------------------

    def get_properties(self):
        
        return
        
# ======================================================================
