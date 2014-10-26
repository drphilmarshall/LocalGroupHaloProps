# ======================================================================

import localgroup

import numpy as np
from sklearn import mixture

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
        self.PDF = -1

        return
        
# ----------------------------------------------------------------------------

    def generate(self,mode="observational",Nsamples=10000):

        self.T.observe_halos(Nsamples=Nsamples)
        self.T.transform_to_M31()
        #dt = np.dtype([('MW_D', 'f8'), ('MW_vr', 'f8'), ('MW_vt', 'f8'), ('M33_D', 'f8'), ('M33_vr', 'f8'), ('M33_vt', 'f8')])
        self.samples = np.transpose(np.array(self.T.get_kinematics()))
        


        # PJM: Might be better to have Triplet.get_kinematics do this 
        # packaging, perhaps... Also, might be better to leave the samples
        # in the Triplet object, and feed them to the GMM...
        
        return
        
# ----------------------------------------------------------------------------

    def approximate(self, mode="GMM"):        
        
        if (mode == "GMM"):
            self.PDF = mixture.GMM()
            self.PDF.fit(self.samples)

        else:
            raise ValueError("Unrecognised approximation mode %s" % mode)

        return

# ======================================================================

if __name__ == '__main__':

    Lhood = Likelihood()
    Lhood.generate()
    print "Sample kinematic parameters: ",Lhood.samples
