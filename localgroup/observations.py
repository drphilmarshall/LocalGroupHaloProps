# ======================================================================

import localgroup

import numpy as np

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
        self.data['LMC'] = dict()
        self.assignment()
        return 

# ----------------------------------------------------------------------------
 
    def assignment(self):
      
        # Milky Way (In heliocentric galactic cartesian):
        self.data['MW']['x'] = [0.00829,0.00016,'Gaussian'] # vdM12, Mpc
        self.data['MW']['y'] = [0.0,0.0,'Delta'] 
        self.data['MW']['z'] = [0.0,0.0,'Delta']
        
        # BUG: need to treat v_rot and v_sun more carefully. See M+C05...
        self.data['MW']['vx'] = [-11.1,1.23,'Gaussian'] # vdM++12 based on Sohn++12 and McMillan12
        self.data['MW']['vy'] = [(-239.0 -12.24),np.sqrt(5.0**2+2.05**2),'Gaussian'] # vdM++12 based on S++12 and McMillan12
        self.data['MW']['vz'] = [-7.25,0.62,'Gaussian'] # vdM++12 based on S++12 and McMillan12

        # M31:
        self.data['M31']['RA'] = [10.6847929,0.08,'Gaussian'] # deg, NED
        self.data['M31']['DEC'] = [41.2690650,0.08,'Gaussian'] # deg, NED
        self.data['M31']['l'] = [121.1744050,0.0,'Delta'] # deg, NED
        self.data['M31']['b'] = [-21.5729360,0.0,'Delta'] # deg, NED
        self.data['M31']['D']  = [0.770,0.040,'Gaussian'] # vdM++12, Mpc

        self.data['M31']['v_west'] = [-125.2, 30.8, 'Gaussian']
        self.data['M31']['v_north'] = [-73.8, 28.4, 'Gaussian']
        
        self.data['M31']['mu_west'] = [-42.2,12.3,'Gaussian'] # Sohn++12, muas/yr, does not contain rotation correction, see vdM++12 for corrections
        self.data['M31']['deltavrot_west'] = [-28.8, 0.0, 'Delta'] # vdM++12, km/s, calculated from numbers in section 2.2
        self.data['M31']['mu_north'] = [-30.9,11.7,'Gaussian'] # S++12, muas/yr, does not contain rotation correction
        self.data['M31']['deltavrot_north'] = [-39.0, 0.0, 'Delta'] # vdM++12, km/s, calculated from numbers in section 2.2

        self.data['M31']['v_r'] = [-301.0,1.0,'Gaussian']  # vdMG08, km/s

        # How were deltavrots derived?
        
        # deltavrot_west  = 4.7404*D*mu_west - v_west 
        # deltavrot_north = 4.7404*D*mu_north - v_north 

        # where v_west, v_north = -125.2 +/- 30.8 , -73.8 +/- 28.4
        #   and D = 0.770, mu_west = -42.2 +/- 12.3, mu_north = -30.9 +/- 11.7
        
        # -> deltavrot_west = -28.8 +/- small: errors dominated by mu uncertainty.
        #    deltavrot_north = -39.0 +/- small
        
        # This is consevrtaive: our uncertainties on v are larger than vdM++12's - because they 
        # use more constraints on the proper motion than Sohn ++2012 even.
        
        # M33:
        self.data['M33']['RA'] = [23.4620417, 0.075,'Gaussian'] # deg, NED
        self.data['M33']['DEC'] = [30.6602222, 0.075,'Gaussian'] # deg, NED
        self.data['M33']['l'] = [133.6101611,0.0,'Delta'] # deg, NED
        self.data['M33']['b'] = [-31.3305875,0.0,'Delta'] # deg, NED
        self.data['M33']['D']  = [0.794,0.023,'Gaussian'] # M++04, vdMG08, Mpc
        self.data['M33']['mu_west'] = [-23.2,6.8,'Gaussian'] # B++05 muas/yr
        self.data['M33']['mu_north'] = [7.5,8.8,'Gaussian'] # B++05 muas/yr
        self.data['M33']['v_r'] = [-180.0,1.0,'Gaussian'] # vdMG08, km/s

        # We may be able to ignore the deltavrot numbers by using the values
        # from Brunthaler++2005, page 4:  mu_west = -23.2 +/- 6.8 (4.7+18.5)
        #                                 mu_north = 7.5 +/- 8.8


        #self.data['LMC']['RA'] = [80.8939, 0.0,'Delta'] # deg, NED
        #self.data['LMC']['DEC'] = [-69.7561, 0.0,'Delta'] # deg, NED

        self.data['LMC']['RA'] = [78.76, 0.52,'Gaussian'] # deg, NED
        self.data['LMC']['DEC'] = [-69.19, 0.25,'Gaussian'] # deg, NED

        self.data['LMC']['l'] = [280.4653,0.0,'Delta'] # deg, NED
        self.data['LMC']['b'] = [-32.8883,0.0,'Delta'] # deg, NED
        self.data['LMC']['D']  = [0.050,1*0.0025,'Gaussian'] # vdM02_LMC Mpc
        self.data['LMC']['mu_west'] = [-1910.0,1*160.0,'Gaussian'] # Kallavayalil13 muas
        self.data['LMC']['mu_north'] = [229.0,1*160.0,'Gaussian'] # Kallavayalil13 muas
        self.data['LMC']['v_r'] = [262.0,3.4,'Gaussian'] # vdMG02_LMC, km/s
        self.data['LMC']['v_west'] = [-399, 40.8, 'Gaussian']
        self.data['LMC']['v_north'] = [80.0, 40.4, 'Gaussian']

        return

# ======================================================================
