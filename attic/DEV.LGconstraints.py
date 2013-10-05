# ======================================================================

import numpy

# ======================================================================

# Construct a dictionary of constraints, for application to the 
# model LG halos. Likelihoods are stored as dictionary entries, 
# {'parameter':[mean,stdev]}

class bundle:

   def __init__(self):
      self.data = dict()
      self.data['MW'] = dict()
      self.data['M31'] = dict()
      self.data['M33'] = dict()
      self.assignment()
      return 
 
# The actual observational data, for named parameters: 
   def assignment(self):

      # # Version 1, radial velocity + distance relative to MW:
      # 
      # # Milky Way:
      # self.data['MW']['Mvir'] = [12.2,0.1,lognormal]
      # 
      # # M31:
      # self.data['M31']['D']  = [0.774,0.040,gaussian] # VG08, Mpc
      # self.data['M31']['vr'] = [-130.0,8.0,gaussian]  # VG08, km/s
      # 
      # # M33:
      # self.data['M33']['D']  = [0.799,0.023,gaussian] # VG08, Mpc
      # self.data['M33']['vr'] = [-190.0,59.0,gaussian] # Brunthaler, km/s
      # 
      # # self.data['M33']['angsep'] = [14,8,0.1,gaussian] # ? Angular separation
      
      # Optional extras - component masses:
      # 
      # Milky Way mass:
      # self.data['MW']['Mvir'] = [12.2,0.1,lognormal] # Busha++2011
      # 
      # M31 mass:
      # self.data['M31']['Mvir'] = [12.4,0.1,lognormal] # ???

      
      # Version 2, radial and tangential velocity, + distance, 
      # relative to M31, all correlated. 
      # Multivariate Gaussian approximation, keep uncertainties in
      # covariance matrix, just assign means here. First need to instantiate 
      # the VanDerMarel class:
      
      VDM = VanDerMarel()

      # Set dictionary of likelihood means:
      self.mean = VDM.mean

      # Set dictionary of indices:
      self.index = VDM.index
      
      # Set up covariance matrix and invert it:
      self.invcovmat = VDM.invcovmat
    
      return


# Given 2 catalogs and a list of parameters to constrain, compute the 
# likelihood of the kinematic data at each point in parameter space.
   def Nd_likelihood(self,catalogs,parameters):

      # Read out columns of catalog into data vectors:
      n = len(catalogs[0]) # No. of sim halos
      # Make parameters the right length:
      np = len(parameters)
      k = len(catalogs)*np # No. of parameters of LG system
      # Predicted data, from catalogs:
      dp = numpy.zeros([n,k])
      # Data, from VDM, and its covariance matrix:
      d = numpy.zeros(k)
      index = [0]*k
      Cinv = numpy.zeros([k,k])
      
      for j,thiscatalog in enumerate(catalogs):
         thisgalaxy = thiscatalog.table_name
         print "Getting observed data for "+thisgalaxy
            
         # Construct the vector of predicted data - sim halo parameters:
         for i,parameter in enumerate(parameters):
            dp[:,i + j*np] = thiscatalog[parameter]
         # This can now be multiplied by an inverse covariance matrix when making
         # naturally weighted residuals.

         # Now construct the corresponding vector of likelihood means:
         for i,parameter in enumerate(parameters):
            d[i + j*np] = self.mean[thisgalaxy][parameter]
            print "Observed value of "+parameter+": ",d[i + j*np]

         # And get the indices that correspond to the parameters:
         for i,parameter in enumerate(parameters):
            index[i + j*np] = self.index[thisgalaxy][parameter]
      
      print "Observed data: ",d
      print "Corresponding indices: ",index
      
      # Finally, pull out the relevant bits of the inv covariance matrix:
      for i in range(k):
         for j in range(k):
            Cinv[i,j] = self.invcovmat[index[i],index[j]]
      print "Inverse covariance matrix for these parameters: ",Cinv
      
      print "Calculating likelihood:"
      
      L = multivariategaussian(d,dp,Cinv)
      
#       x = catalog[parameter]
#       x0 = d[0]
#       sigma = d[1]
#       func = d[2]
#       L = apply(func,(x,x0,sigma))
      
      return L

# ======================================================================
# Global functions needed by the bundle class - unit maximum PDFs:

def gaussian(x,x0,sigma):
   return numpy.exp(-(x-x0)*(x-x0)/(2.0*sigma*sigma))

def halfgaussian(x,x0,sigma):
# Sigma is x within which 68% of the probability is included: 
# \int_0.0^sigma halfgaussian dx = 0.68* \int_0.0^infty halfgaussian dx
# ie its the same sigma, by symmetry!
   return numpy.exp(-(x-x0)*(x-x0)/(2.0*sigma*sigma))

def lognormal(x,x0,sigma):
   y = numpy.log10(x)
   return gaussian(y,x0,sigma)
   
def multivariategaussian(x,x0,Cinv):
   assert Cinv.shape[0] == Cinv.shape[1]
   print "mvg: x = ",x,x.shape
   print "mvg: x0 = ",x0,x0.shape
   print "mvg: Cinv = ",Cinv,Cinv.shape
# This needs to be a matrix multiply... broadcasting error:
# ERROR: ValueError: operands could not be broadcast together with shapes (2583,3) (3,3)  [LGconstraints]   
   chisq = (x-x0)*Cinv*(x-x0)
# chisq needs to be a vector of length np...   
   print "mvg: chisq = ",chisq,chisq.shape
   return numpy.exp(-chisq/2.0)

# ======================================================================
# Class to handle the galactocentric postion and velocity data of 
# VdM12. Given mean and covmat for x,y,z,vx,vy,vz, compute mean and covmat
# for D,vr,vt by Monte Carlo. Save estimate in pickle.

class VanDerMarel:

      def __init__(self):
         self.set_mean()
         self.set_index()
         self.set_invcovmat()
         return 
 
      def set_index(self):
         self.index = dict()
         self.index['MW'] = dict()
         self.index['MW']['D']  = 0
         self.index['MW']['vr'] = 1
         self.index['MW']['vt'] = 2
         self.index['M33'] = dict()
         self.index['M33']['D']  = 3
         self.index['M33']['vr'] = 4
         self.index['M33']['vt'] = 5
         self.nd = 6
         return

      def set_mean(self):
         self.mean = dict()
         self.mean['MW'] = dict()
         self.mean['MW']['D']  =  0.774 # Mpc
         self.mean['MW']['vr'] = -130.0 # km/s
         self.mean['MW']['vt'] =    0.0 # km/s
         self.mean['M33'] = dict()
         self.mean['M33']['D']  = 0.300 # Mpc
         self.mean['M33']['vr'] = -50.0 # km/s
         self.mean['M33']['vt'] =   0.0 # km/s
         return
         
      def set_invcovmat(self):
         self.covmat = numpy.zeros([self.nd,self.nd])
         self.covmat[0,0] =  0.04 # Mpc
         self.covmat[1,1] =   8.0 # km/s
         self.covmat[2,2] =  50.0 # km/s
         self.covmat[3,3] =  0.04 # Mpc
         self.covmat[4,4] =  16.0 # km/s
         self.covmat[5,5] = 100.0 # km/s
         self.invcovmat = numpy.linalg.inv(self.covmat)
         return

# ======================================================================
# \hline
# $\rx$ / kpc           & $-378.9 \pm 30.?$ & $-476.1 \pm 30.?$ & vdM12   \\
# $\ry$ / kpc           & $ 612.7 \pm 30.?$ & $ 491.1 \pm 30.?$ & vdM12   \\
# $\rz$ / kpc           & $-283.1 \pm 30.?$ & $-412.9 \pm 30.?$ & vdM12   \\
# $\vx$ / km s$^{-1}$   & $  66.1 \pm 26.7$ & $  43.1 \pm 21.3$ & vdM12   \\
# $\vx$ / km s$^{-1}$   & $ -76.3 \pm 19.0$ & $ 101.3 \pm 23.5$ & vdM12   \\
# $\vx$ / km s$^{-1}$   & $  45.1 \pm 26.5$ & $ 138.8 \pm 28.1$ & vdM12   \\
# \hline
#                       & MW                & M33               & \\
# \hline
# $\distance$ / kpc     & $ xxx.x \pm xx.x$ & $ xxx.x \pm xx.x$ & \\
# $\vrad$ / km s$^{-1}$ & $  xx.x \pm xx.x$ & $  xx.x \pm xx.x$ & \\
# $\vtan$ / km s$^{-1}$ & $  xx.x \pm xx.x$ & $  xx.x \pm xx.x$ & \\
# ======================================================================

# Testing!

if __name__ == '__main__':
  
# Make an instance of a VanDerMarel2012 Gaussian:
  VDM = VanDerMarel()

# Draw samples in {{D,vr,vt}_MW,{D,vr,vt}_M33} from this, and write out 
# for plotting!
  ns = 1000
  samples = numpy.random.multivariate_normal(VDM.mean,VDM.covmat,ns)
  print samples.shape, samples


# ======================================================================
