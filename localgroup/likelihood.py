# ======================================================================

import localgroup

import numpy as np
from sklearn import mixture
from sklearn.grid_search import GridSearchCV
from scipy import linalg
import matplotlib as mpl
import triangle
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
        generate           - draw sample vector from observations' distribuions
        approximate        - compute KNN/GMM/etc estimate of PDF
        evaluate           - compute value of PDF at given vector
        plot_samples       - make triangle plot for observed data in M31 ref frame
        set_PDF            - set the L.PDF field
        test_gauss         - plots and calculates score vs ngauss components
        preprocess_samples - zero the mean of the samples and scale by standard deviation

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
        self.PDF = None

        return
        
# ----------------------------------------------------------------------------

    def generate(self,mode="observational",Nsamples=10000):

        self.T.observe_halos(Nsamples=Nsamples)
        self.T.transform_to_M31()
        #dt = np.dtype([('MW_D', 'f8'), ('MW_vr', 'f8'), ('MW_vt', 'f8'), ('M33_D', 'f8'), ('M33_vr', 'f8'), ('M33_vt', 'f8')])
        self.samples = np.transpose(np.array(self.T.get_kinematics()))
        self.samples_means = np.array([np.mean(self.samples[:,i]) for i in range(self.samples.shape[1])])
        self.samples_stds = np.array([np.std(self.samples[:,i]) for i in range(self.samples.shape[1])])

        self.preprocess_samples()


        # PJM: Might be better to have Triplet.get_kinematics do this 
        # packaging, perhaps... Also, might be better to leave the samples
        # in the Triplet object, and feed them to the GMM...
        
        return
        
# ----------------------------------------------------------------------------

    def approximate(self, mode="GMM", cv=False):        
        
        if (mode == "GMM"):
            if (cv):
                def bic_scorefn(GMM, X): return GMM.bic(X)
                score_dict = self.test_GaussMM(num_folds=10, score_fns=[bic_scorefn], maxMM=15)
                bic = score_dict['bic_scorefn']

                min_bic = np.min(bic)
                min_bic_components = 1+np.array([i for i, x in enumerate(bic) if x == min_bic])
                print "Minimum BIC: ", min_bic
                print "Minimum BIC number of components: ", min_bic_components
                self.set_PDF(mixture.GMM(min_bic_components[0], covariance_type='full'))
            self.PDF.fit(self.samples)
        else:
            raise ValueError("Unrecognised approximation mode %s" % mode)

        return
# ======================================================================

    def evaluate(self, points, mode="GMM"):
        if (mode == "GMM"):
            logprobs, resps = self.PDF.score_samples(points)
        else:
            raise ValueError("Unrecognised approximation mode %s" % mode )

        return logprobs, resps

# ======================================================================

    def plot_samples(self, ngauss, overlay=False):
        try:
            figure = triangle.corner(self.samples, labels=["MW_D", "MW_vr", "MW_vt", "M33_D", "M33_vr", "M33_vt"], quantiles=[0.16,0.5,0.84], show_titles=True, title_args={"fontsize": 12})
        except AttributeError:
            raise AttributeError("L.generate has not been run.")
        figure.gca().annotate("MW and M33 Observational Data Distributions (M31 centric)", xy=(0.5, 1.0), xycoords="figure fraction", xytext=(0, -5), textcoords="offset points", ha="center", va="top")
        
        if overlay:
            self.gaussianOverlay(figure, ngauss)
        figure.savefig("L_samples_tri.png")
        return figure

# ======================================================================

    def set_PDF(self, pdf):
        self.PDF = pdf
        return

# ======================================================================

    def test_GaussMM(self, num_folds, score_fns, maxMM=10):
        aic_scores = []
        bic_scores = []
        scores_dict = {}
        params = {'n_components':np.arange(1,maxMM+1), 'covariance_type':['full']}
        gmm = mixture.GMM()
        for score_fn in score_fns:
            grid = GridSearchCV(estimator=gmm, param_grid=params, cv=num_folds, scoring=score_fn)
            grid.fit(self.samples)
            scores_dict[score_fn.func_name] = np.array(grid.grid_scores_)[:,1]

        return scores_dict

# ======================================================================

    def preprocess_samples(self):

        self.samples = (self.samples - self.samples_means)/self.samples_stds

        return

# ======================================================================


    def gaussianOverlay(self, figure, ngauss):
        n_gaussians = ngauss
        if ngauss > 5: raise AttributeError("Only 5 colors can be shown.")
        colors = ['g', 'r', 'y', 'b', 'c']
        transparency = 0.5
        model = mixture.GMM(n_gaussians, covariance_type='full')
        axes = np.reshape(figure.axes, (6,6))
        for i in range(6):
            for j in range(6):
                if j < i:
                    model.fit(self.samples)
                    subplot = axes[i,j]
                    for gauss_num in range(n_gaussians):
                        mean = [model.means_[gauss_num][j], model.means_[gauss_num][i]]
                        covar = model.covars_[gauss_num]
                        covar = [[covar[j,j], covar[j,i]], [covar[i,j], covar[i,i]]]
                        color = colors[gauss_num]
                        v, w = linalg.eigh(covar)
                        u = w[0]/linalg.norm(w[0])
                        angle = np.arctan(u[1]/u[0])
                        angle = 180 * angle / np.pi
                        ell = mpl.patches.Ellipse(mean, 2*np.sqrt(v[0]), 2*np.sqrt(v[1]), 180 + angle, color=color)
                        ell.set_clip_box(subplot.bbox)
                        ell.set_alpha(transparency)
                        subplot.add_artist(ell)

        return




# ======================================================================


if __name__ == '__main__':

    Lhood = Likelihood()
    Lhood.generate()
    print "Sample kinematic parameters: ",Lhood.samples
