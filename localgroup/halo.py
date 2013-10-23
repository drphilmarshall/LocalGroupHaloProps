# ======================================================================

import localgroup

import numpy

deg2rad = numpy.pi/180.0
maspyrMpc2kmps = 4.7404

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

    def sample_from(self,obs,Nsamples):
        
        # import astropy.units as u
        # import astropy.coordinates as coords
        
        self.RA              = localgroup.draw(obs['RA'],Nsamples)
        self.DEC             = localgroup.draw(obs['DEC'],Nsamples)
        self.D               = localgroup.draw(obs['D'],Nsamples)
        self.mu_west         = localgroup.draw(obs['mu_west'],Nsamples)
        self.deltavrot_west  = localgroup.draw(obs['deltavrot_west'],Nsamples)
        self.mu_north        = localgroup.draw(obs['mu_north'],Nsamples)
        self.deltavrot_north = localgroup.draw(obs['deltavrot_north'],Nsamples)
        self.vr              = localgroup.draw(obs['vr'],Nsamples)
        
        # Equatorial tangential motions in km/s:
        self.v_west  = maspyrMpc2kmps*self.D*self.mu_west - self.deltavrot_west
        self.v_north = maspyrMpc2kmps*self.D*self.mu_north - self.deltavrot_north

        # # Heliocentric cartesian positions:
        # self.x = numpy.zeros(Nsamples)
        # self.y = numpy.zeros(Nsamples)
        # self.z = numpy.zeros(Nsamples)
        # for i,(ra,dec,D) in enumerate(zip(self.RA, self.DEC, self.D)):
        #     icrs = coords.ICRSCoordinates(ra, dec, unit=(u.deg, u.deg),
        #                       distance=coords.Distance(D,u.Mpc))
        #     self.x[i] = icrs.x
        #     self.y[i] = icrs.y
        #     self.z[i] = icrs.z
        # This loop is slow - better to do transformations by hand!

        # Heliocentric cartesian positions (Mpc):
        delta = self.DEC*deg2rad
        alpha = self.RA*deg2rad
        self.x = self.D*numpy.cos(delta)*numpy.cos(alpha)
        self.y = self.D*numpy.cos(delta)*numpy.sin(alpha)
        self.z = self.D*numpy.sin(delta)

        # # Heliocentric cartesian velocities (km/s):
        # self.vx = (self.x/self.D)*self.vr - (self.y/self.D)*self.v_west - (self.z/self.D)*numpy.cos(alpha)*self.v_north
        # self.vy = (self.y/self.D)*self.vr + (self.x/self.D)*self.v_west - (self.z/self.D)*numpy.sin(alpha)*self.v_north
        # self.vz = (self.z/self.D)*self.vr + numpy.cos(delta)*self.v_north
         
        # Heliocentric cartesian velocities (km/s):
        self.vx = (self.x/self.D)*self.vr - (self.z/self.D)*self.v_north*numpy.cos(alpha) - self.v_west *numpy.sin(alpha) 
        self.vy = (self.y/self.D)*self.vr - (self.z/self.D)*self.v_north*numpy.sin(alpha) + self.v_west *numpy.cos(alpha) 
        self.vz = (self.z/self.D)*self.vr + numpy.cos(delta)*self.v_north

        print "Position: ",self.x[0],self.y[0],self.z[0]
        print "Velocity (xyz): ",self.vx[0],self.vy[0],self.vz[0]
        print "Velocity (rWN): ",self.vr[0],self.v_west[0],self.v_north[0]

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
        
        if newsystem == 'Heliocentric Cartesians':
            self.coordinate_system == newsystem
            pass
        
        return
        
# ======================================================================
