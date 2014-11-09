# ======================================================================

import localgroup
import numpy as np
from scipy import linalg

# ======================================================================

class GaussMM(object):
    """
    NAME
        GaussMM

    PURPOSE
        Create a Gaussian Mixture Model object.

    COMMENTS
        This class allows you to create a Gaussian Mixture Model with n components.
    INITIALISATION

    
    METHODS
        fit - Trains a GMM with n components using the EM algorithm.

    BUGS

    AUTHORS
      This file is part of the LocalGroupHaloProps project, 
      distributed under the GPL v2, 
      by Marc Williamson Phil Marshall (KIPAC). 
      Please cite: Williamson et al in preparation.

    HISTORY
      2014-11-08  started Williamson and Marshall (KIPAC)
    """
# ======================================================================

    def __init__(self, ncomponents):
        
        self.ngauss = ncomponents
        self.means = None
        self.covars = None
        self.phi = 1.0/ncomponents*np.ones(ncomponents)

        return
        
# =====================================================================

    def fit(self, train_data, EM_cut_off=50):
        num_EM_steps = 0
        continue_training = True
        num_rows, num_feat = train_data.shape
        self.means = np.ones((self.ngauss, num_feat))
        self.covars = np.ones((self.ngauss, num_feat, num_feat))
        while(continue_training and num_EM_steps < EM_cut_off):
            weights = self.setweights(train_data, num_rows, num_feat)
            for j in range(self.ngauss):
                self.phi[j] = 1.0/num_rows*np.sum(weights[:,j])
                self.means[j] = 1.0/np.sum(weights[:,j])*weights[:,j].dot(train_data)
                temp_cov = np.zeros((num_feat, num_feat))
                for i in range(num_rows):
                    temp_cov = temp_cov + np.array(weights[i,j]*np.matrix(train_data[i]-self.means[j]).transpose().dot(np.matrix(train_data[i]-self.means[j])))
                self.covars[j] = temp_cov
            num_EM_steps = num_EM_steps + 1
        return

# ====================================================================

    def setweights(self, train_data, num_rows, num_feat):

         weights = np.ndarray((num_rows, self.ngauss))
         for i in range(num_rows):
             for j in range(self.ngauss):
                 weights[i,j] = 1.0/( (2*np.pi)**(num_feat/2.0)*np.sqrt(linalg.det(self.covars[j])) )
                                  *np.exp(-1.0/2.0*train_data[i].dot(linalg.pinv(self.covars[j])).dot(train_data[i]))

         return weights

















if __name__ == '__main__':

    Lhood = Likelihood()
    Lhood.generate()
    print "Sample kinematic parameters: ",Lhood.samples
