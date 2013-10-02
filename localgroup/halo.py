# ======================================================================

import localgroup

import numpy

# ======================================================================

class Halo(object):
    """
    NAME
        Halo

    PURPOSE
        Define a dark matter halo (or set of halos), that has either:
         a) observables, like sky position, proper motion etc, or
         b) cartesian phase space coordinates, {x,y,z,vx,vy,vz}
        and potentially
         c) properties, like Mass, Angular Momentum, etc,
        and has the means to translate between its coordinate systems, 
        transform observables to coordinates, and compute its likelihood
        of being a local group analog.

    COMMENTS
        This class is used in two places: derivation of the local group
        analog likelihood function, and then inference of local group
        halo properties using that likelihood.

    INITIALISATION

    
    METHODS

        read
        sample_from
        translate_to
        transform_to

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

    def read(self,halofile):
        
        self.properties = astropy.Table(halofile)
        
        self.x = self.properties['x']
        self.y = self.properties['y']
        self.z = self.properties['z']
        self.vx = self.properties['vx']
        self.vy = self.properties['vy']
        self.vz = self.properties['vz']

        return
        
# ----------------------------------------------------------------------------

    def sample_from(self,obs,Nsamples):
        
        self.ra = draw(obs['RA'],Nsamples)
        # etc
        
        return
        
# ----------------------------------------------------------------------------

    def translate_to(self,other):
        
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        self.vx -= other.vx
        self.vy -= other.vy
        self.vz -= other.vz
        self.frame = other.name
        
        return
        
# ----------------------------------------------------------------------------

    def transform_to(self,newsystem):
        
        if newsystem == 'Heliocentric Cartesians' then:
            self.coordinate_system == newsystem
            pass
        
        return
        
# ======================================================================
