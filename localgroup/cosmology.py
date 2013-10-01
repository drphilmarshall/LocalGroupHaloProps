# ======================================================================

# Useful cosmology functions.

import numpy

# Default parameters:

Om0 = 0.3
OL0 = 0.7
h = 0.7

G = 4.30e-9 # km^2 s^-2 Mpc Msun^-1

# ======================================================================
# Critical density at redshift z:

def rhocrit(z):
   return 3.0*Hsq(z)/(8.0*numpy.pi*G)  # Msun / Mpc^-3

# ----------------------------------------------------------------------
# Hubble function in non-flat, Lambda CDM universe:      

def Hsq(z):
   return Esq(z)*(100.0*h)**2

# ----------------------------------------------------------------------
# Dimensionless Hubble function:      

def Esq(z):
   return (Om0*(1+z)**3 + OL0 + (1-Om0-OL0)*(1+z)**2)

# ----------------------------------------------------------------------
# Matter density parameter at redshift z:      

def Om(z):
   return Om0*(1+z)**3/Esq(z)

# ----------------------------------------------------------------------
# Virial overdensity at redshift z:      

def Deltav(z):
   return Deltac(z)/Om(z)

# ----------------------------------------------------------------------
# Critical overdensity at redshift z:      

def Deltac(z):
   x = Om(z) - 1.0
   return 18*numpy.pi**2 + 82*x + 39*x**2

# ======================================================================     

# Testing:

if __name__ == '__main__':

   z = numpy.array([0.0,0.0,0.0])
   
   print "z = ",z
   print "rhocrit(z) = ",rhocrit(z)
   print "Hsq(z) = ",Hsq(z)
   print "Om(z) = ",Om(z)
