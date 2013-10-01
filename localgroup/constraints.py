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


      # Version 1, radial velocity + distance relative to MW:
      
      # Milky Way:
      self.data['MW']['Mvir'] = [12.2,0.1,lognormal]
      
      # M31:
      self.data['M31']['D']  = [0.774,0.040,gaussian] # VG08, Mpc
      self.data['M31']['vr'] = [-130.0,8.0,gaussian]  # VG08, km/s
      
      # M33:
      self.data['M33']['D']  = [0.799,0.023,gaussian] # VG08, Mpc
      self.data['M33']['vr'] = [-190.0,59.0,gaussian] # Brunthaler, km/s
      
      # self.data['M33']['angsep'] = [14,8,0.1,gaussian] # ? Angular separation
      
      
#       # Version 2, radial velocity + distance, relative to M31:
# 
#       # MW, relative to M31:
#       # NEED TO UPDATE TO VdM++2012!
#       self.data['MW']['D']  = [0.774,0.040,gaussian]  # VG08, Mpc
#       self.data['MW']['vr'] = [-130.0,8.0,gaussian]   # VG08, km/s
#       # self.data['MW']['vt'] = [0.0,4.0,halfgaussian]  # VG12, km/s
# 
#       # M33, relative to M31:
#       # NEED TO COMPUTE FROM VdM++2012!
#       self.data['M33']['D']  = [0.799,0.023,gaussian] # VG08, Mpc
#       self.data['M33']['vr'] = [-190.0,59.0,gaussian] # Brunthaler, km/s
#       # self.data['M33']['vt'] = [0.0,4.0,halfgaussian]  # To be calculated!!
# 
#       # Optional extras - component masses:
# 
#       # Milky Way mass:
#       self.data['MW']['Mvir'] = [12.2,0.1,lognormal] # Busha++2011
#       
#       # M31 mass:
#       self.data['M31']['Mvir'] = [12.4,0.1,lognormal] # ???

      return


   def likelihood(self,catalog,parameter):

      thisgalaxy = catalog.table_name
      d = self.data[thisgalaxy][parameter]
      print "Calculating likelihoods: galaxy "+thisgalaxy+", data = ",d
      
      x = catalog[parameter]
      x0 = d[0]
      sigma = d[1]
      func = d[2]
      L = apply(func,(x,x0,sigma))
      
      return L

# ======================================================================
# Global functions needed by the class - unit maximum PDFs:

def gaussian(x,x0,sigma):
   return numpy.exp(-(x-x0)*(x-x0)/(2.0*sigma*sigma))

def lognormal(x,x0,sigma):
   y = numpy.log10(x)
   return gaussian(y,x0,sigma)   

# ======================================================================
