# ======================================================================

import localgroup

import numpy

# ======================================================================

class Likelihood(object):
    """
    NAME
        Likelihood

    PURPOSE
        Compute the likelihood of a given halo having the observed properties.

    COMMENTS
        This class allows you to generate the approximation used, as well as evalu        ate the likelihood. 

    INITIALISATION

    
    METHODS
        generate         - draw sample vector from observations' distribuions
        approximate      - compute KNN/GMM/etc estimate of PDF
        evaluate         - compute value of PDF at given vector
        
        NB. "vector" refes to a position in 6D MW-M31 D,vr,vt space

    BUGS

    AUTHORS
      This file is part of the LocalGroupHaloProps project, 
      distributed under the GPL v2, 
      by Marc Williamson Phil Marshall (KIPAC). 
      Please cite: Williamson et al in preparation.

    HISTORY
      2014-09-23  started Williamson and Marshall (KIPAC)
    """
# ======================================================================

    def __init__(self):
        
        self.T = localgroup.Triplet()
    
        return
        
# ----------------------------------------------------------------------------

    def generate(self,mode="observational"):

        self.T.observe_halos()
        
                
        return
        
# ----------------------------------------------------------------------------
