# ======================================================================

import localgroup

import numpy

deg2rad = numpy.pi/180.0
muaspyrMpc2kmps = 4.7404

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

    def __init__(self,name):
        
        self.name = name
        
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
# Draw from sampling distributions associated with measurements, mostly
# in equatorial coordinates (for external galaxies) and convert to 
# heliocentric galactic cartesian coordinates.

    def sample_from(self,obs,Nsamples):
                
        if self.name == 'MW':
            
            # Sample from Cartesian coordinates directly:
            self.x               = localgroup.draw(obs['x'],Nsamples)
            self.y               = localgroup.draw(obs['y'],Nsamples)
            self.z               = localgroup.draw(obs['z'],Nsamples)
            self.vx              = localgroup.draw(obs['vx'],Nsamples)
            self.vy              = localgroup.draw(obs['vy'],Nsamples)
            self.vz              = localgroup.draw(obs['vz'],Nsamples)
                       
        else:
         
            # Get observed quantities and transform:
            self.RA              = localgroup.draw(obs['RA'],Nsamples)
            self.DEC             = localgroup.draw(obs['DEC'],Nsamples)
            self.l               = localgroup.draw(obs['l'],Nsamples)
            self.b               = localgroup.draw(obs['b'],Nsamples)
            self.D               = localgroup.draw(obs['D'],Nsamples)
            self.mu_west         = localgroup.draw(obs['mu_west'],Nsamples)
            self.deltavrot_west  = localgroup.draw(obs['deltavrot_west'],Nsamples)
            self.mu_north        = localgroup.draw(obs['mu_north'],Nsamples)
            self.deltavrot_north = localgroup.draw(obs['deltavrot_north'],Nsamples)
            self.v_r             = localgroup.draw(obs['v_r'],Nsamples)

            # Equatorial tangential motions in km/s:
            self.v_west  = muaspyrMpc2kmps*self.D*self.mu_west - self.deltavrot_west
            self.v_north = muaspyrMpc2kmps*self.D*self.mu_north - self.deltavrot_north

            # Convert proper motions to galactic coordinates:
            self.v_l,self.v_b = localgroup.equatorial_to_galactic_proper_motion(self.v_west,self.v_north,self.RA,self.DEC)

            # Convert from spherical to cartesian coordinates:
            self.x,self.y,self.z,self.vx,self.vy,self.vz = \
                localgroup.spherical_to_cartesian(self.l,self.b,self.D,self.v_l,self.v_b,self.v_r)
                
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
                
# ======================================================================
