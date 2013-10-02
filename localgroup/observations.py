# ======================================================================

import localgroup

import numpy

# ======================================================================

class Observations(object):
    """
    NAME
        Observations

    PURPOSE
        Encode all the observational information we have about the 
        kinematics of the 3 largest local group galaxies.

    COMMENTS

    INITIALISATION

    METHODS

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
        self.data = dict()
        self.data['MW'] = dict()
        self.data['M31'] = dict()
        self.data['M33'] = dict()
        self.assignment()
        return 

# ----------------------------------------------------------------------------
 
    def assignment(self):
      
        # Milky Way:
        self.data['MW']['RA'] = [value,err,func]

        # M31:
        self.data['M31']['D']  = [0.774,0.040,'Gaussian'] # VG08, Mpc
        self.data['M31']['vr'] = [-130.0,8.0,'Gaussian']  # VG08, km/s

        # M33:
        self.data['M33']['D']  = [0.799,0.023,'Gaussian'] # VG08, Mpc
        self.data['M33']['vr'] = [-190.0,59.0,'Gaussian'] # Brunthaler, km/s

        return

# ======================================================================
