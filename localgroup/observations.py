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
      
#         # Milky Way:
#         self.data['MW']['RA'] = [value,err,'Gaussian']
#         self.data['MW']['DEC'] = [value,err,'Gaussian']
#         self.data['MW']['D'] = [0.00829,0.00016,'Gaussian'] #vdM12, Mpc
#         self.data['MW']['mu_west'] = [value,err,'Gaussian'] 
#         self.data['MW']['mu_north'] = [value,err,'Gaussian']
#         self.data['MW']['vr'] = [value,err,'Gaussian']

        # M31:
        self.data['M31']['RA'] = [10.6847929,0.0,'Delta'] #deg, NED
        self.data['M31']['DEC'] = [41.2690650,0.0,'Delta'] #deg, NED
        self.data['M31']['D']  = [0.774,0.040,'Gaussian'] # vdMG08, Mpc
        self.data['M31']['mu_west'] = [-0.0422,0.0123,'Gaussian'] #S++12, mas/yr, does not contain rotation correction, see vdM++12 for corrections
        self.data['M31']['deltavrot_west'] = [8.7, 13.9, 'Gaussian'] #vdM++12, km/s, calculated from numbers in section 2.2
        self.data['M31']['mu_north'] = [-0.0309,0.0117,'Gaussian'] #S++12, mas/yr, does not contain rotation correction
        self.data['M31']['deltavrot_north'] = [4.3, 13.6, 'Gaussian'] #vdM++12, km/s, calculated from numbers in section 2.2
        self.data['M31']['vr'] = [-301.0,1.0,'Gaussian']  # vdMG08, km/s


        # M33:
        self.data['M33']['RA'] = [23.4620417, 0.0,'Delta'] #deg, NED
        self.data['M33']['DEC'] = [30.6602222, 0.0,'Delta'] #deg, NED
        self.data['M33']['D']  = [0.794,0.023,'Gaussian'] # M++04, vdMG08, Mpc
        self.data['M33']['mu_west'] = [-4.7,3.2,'Gaussian'] #vDMG08, B++05, B++07 mas/yr
        self.data['M33']['deltavrot_west'] = [70.0, 23.0, 'Gaussian'] #vdMG08, km/s, derived from B++05
        self.data['M33']['mu_north'] = [-14.1,6.4,'Gaussian'] #vDMG08, B++05, B++07 mas/yr
        self.data['M33']['deltavrot_north'] = [-81.0, 23.0, 'Gaussian'] #vdMG08, km/s, derived from B++05
        self.data['M33']['vr'] = [-180.0,1.0,'Gaussian'] # vdMG08, km/s

        return

# ======================================================================
