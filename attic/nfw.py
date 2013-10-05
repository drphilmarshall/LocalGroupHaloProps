# ======================================================================

import numpy
import cosmology

vb = False

# ======================================================================

# Compute various quantities for the NFW halo model.

class halo:

   def __init__(self,M,z,c=None,kind='virial'):
      
      self.z = z
      
      if kind == 'virial':
      # Compute quantities at R200 from virial ones:
         self.Mv = M
         if c != None: 
           self.cv = c
         else:
           self.mcrelation(kind)  
         
      elif kind == '200c':
         self.M200 = M
         if c != None: 
           self.c200 = c
         else:
           self.mcrelation(kind)

      else:
         raise "Unrecognised mass kind "+kind

      # Compute R, rs, rhos for completeness:
      self.otherpars(kind)

      return

   # -------------------------------------------------------------------   
   #  Have M of some kind, get corresponding c.
   
   def mcrelation(self,kind):
   
      if kind == 'virial':
      # Transformed from Neto et al 2007:
         self.cv = 5.54416294*(cosmology.h*self.Mv/1e14)**(-0.11) # placeholder
         
      elif kind == '200c':
      # Neto et al 2007, Millenium simulations, all halos, 12 < log M < 15:
         self.c200 = 4.67*(cosmology.h*self.M200/1e14)**(-0.11)
      
      else:
         raise "Unrecognised mass kind "+kind

      return
      
   # -------------------------------------------------------------------   
   #  Have M and c of some kind, get R, rs and rhos.
   
   def otherpars(self,kind):
   
      if kind == 'virial':
         self.Rv = ((3*self.Mv)/(4*numpy.pi*cosmology.Deltac(self.z)*cosmology.rhocrit(self.z)))**(1.0/3.0)
         self.rs = self.Rv/self.cv
         self.rhos = self.Mv/(4*numpy.pi*self.Rv**3*hk_f(1/self.cv))
         
      elif kind == '200c':
         self.R200 = ((3*self.M200)/(4*numpy.pi*200*cosmology.rhocrit(self.z)))**(1.0/3.0)
         self.rs = self.R200/self.c200
         self.rhos = self.M200/(4*numpy.pi*self.R200**3*hk_f(1/self.c200))
      
   # -------------------------------------------------------------------   
   #  Return value of density (in Msun Mpc^-3) at radius r.
   #  Requires rs and rhos.
   
   def density(self,r):
      x = r / self.rs
      return self.rhos/(x*(1.0+x)**2)
      
   # -------------------------------------------------------------------   
   #  Return value of surface density (in Msun Mpc^-2) at projected radius
   #  R. Requires rs and rhos.
   
   def surface_density(self,R):
      x = R / self.rs
      return (self.rhos*self.rs)*sf(x)
      
   # -------------------------------------------------------------------   
   #  Return value of dimensioned shear (in Msun Mpc^-2) at projected radius
   #  R. Requires rs and rhos.
   
   def shear(self,R):
      x = R / self.rs
      return (self.rhos*self.rs)*sg(x)
      
   # -------------------------------------------------------------------
   # Have Mv and cv, want M200 and c200. Follow Hu & Kravtsov:
  
   def massconvert(self,kind):
      
      if kind == 'virial':
      
        fv = hk_f(1/self.cv)
        f200 = 200.0*fv/cosmology.Deltac(self.z)
        x = hk_x(f200) # = rs/R200, note this is c200, upside-down
        
        self.c200 = 1.0/x        
        Rratio = self.c200/self.cv
        self.M200 = self.Mv*(f200/fv)*Rratio**3
        self.R200 = Rratio*self.Rv
      
      elif kind == '200c':
      
        f200 = hk_f(1/self.c200)
        fv = f200*cosmology.Deltac(self.z)/200.0
        x = hk_x(fv)
        
        self.cv = 1.0/x
        Rratio = self.cv/self.c200
        self.Mv = self.M200*(fv/f200)*Rratio**3
        self.Rv = Rratio*self.R200

      else:
      
         raise "Unrecognised mass kind "+kind
      
      return
      
      
# ======================================================================     
# Global functions needed by massconvert:

def hk_f(x):
   return x*x*x*(numpy.log(1.0 + 1.0/x) - 1.0/(1.0+x))
   
hk_a = [0.0,0.5116,-0.4283,-3.13e-3,-3.52e-5]
def hk_x(f):
   logf = numpy.log(f)
   twop = 2*(hk_a[2] + hk_a[3]*logf + hk_a[4]*logf*logf)
   return 1.0/numpy.sqrt(hk_a[1]*f**twop + 0.75**2) + 2*f

# ======================================================================     
# Global functions needed by surface density etc:

def sf(x):
   if x>1:
      z = (1-(2./(x**2-1)**.5)*numpy.arctan(((x-1.)/(x+1))**.5))/(x**2-1.)
   elif x==1:
      z = 1.0/3.0
   else:
      z = (1.-(2./(1-x**2)**.5)*numpy.arctanh(((1.-x)/(x+1))**.5))/(x**2-1)
   return 2*z

def sg(x):
   if x>1:
      y = (((x-1)/(x+1))**0.5)
      z = (8.0*numpy.arctan(y)  / (x**2*(x**2-1.0)**0.5)) +\
          (4.0/x**2)*numpy.log(x/2.0) - \
           2.0/(x**2-1.0) +\
           4.0*numpy.arctan(y)  / (((x**2)-1.0)**(3.0/2.0))
   elif x==1:
      z =(10.0/3.0 + 4*numpy.log(0.5))
   else:
      y = (((1-x)/(x+1))**0.5)
      z = (8.0*numpy.arctanh(y) / (x**2*(1.0-x**2)**0.5)) +\
          (4.0/x**2)*numpy.log(x/2.0) - \
           2.0/(x**2-1.0) +\
           4.0*numpy.arctanh(y) / ((x**2-1.0)*(1.0-x**2)**(1.0/2.0))
   return z
   
# ======================================================================     

# Testing:

if __name__ == '__main__':

#    M = numpy.array([1e12,1e13,1e14])
#    # z = numpy.array([0.0,0.2,0.4])
#    z = 0.0
   M = 2.3e13
   z = 0.1
   vb = True
   
#    h = halo(M,z,kind='virial')   
# 
#    print "          Halo redshift z = ",z
#    print "   rho_crit(z)/MsunMpc^-3 = ",cosmology.rhocrit(z)
#    print "Virial overdensity Deltav = ",cosmology.Deltav(z)
#    print "Virial overdensity Deltac = ",cosmology.Deltac(z)
#    print " "
#    print "      Virial Mass Mv/Msun = ",h.Mv
#    print "  Virial concentration cv = ",h.cv
#    print "     Virial radius Rv/Mpc = ",h.Rv
#    print " "
#    print "      Scale radius rs/Mpc = ",h.rs
#    rhosv = h.rhos
#    print "    Scale rhos/MsunMpc^-3 = ",rhosv
#    print " "
#    h.massconvert('virial')
#    print "                M200/Msun = ",h.M200
#    print "       concentration c200 = ",h.c200
#    print "     Halo radius R200/Mpc = ",h.R200
#    print " "
#    h.otherpars('200c')
#    print "      Scale radius rs/Mpc = ",h.rs
#    print "    Scale rhos/MsunMpc^-3 = ",h.rhos
#    print "     Halo radius R200/Mpc = ",h.R200
#    print " "
#    print "Mass difference (%) = ",100*(h.Mv-h.M200)/h.Mv
#    print "rhos difference (%) = ",100*(rhosv-h.rhos)/rhosv
   
   
   # Invert!
#    i = halo(h.M200,z,kind='200c')   
   i = halo(M,z,kind='200c')   

   print "          Halo redshift z = ",z
   print "   rho_crit(z)/MsunMpc^-3 = ",cosmology.rhocrit(z)
   print "Virial overdensity Deltav = ",cosmology.Deltav(z)
   print "Virial overdensity Deltac = ",cosmology.Deltac(z)
   print " "
   print "           Mass M200/Msun = ",i.M200
   print "       concentration c200 = ",i.c200
   print "     Halo radius R200/Mpc = ",i.R200
   print " "
   print "      Scale radius rs/Mpc = ",i.rs
   rhos = i.rhos
   print "    Scale rhos/MsunMpc^-3 = ",rhos
   print " "

#    i.massconvert('200c')
#    print "      Virial Mass Mv/Msun = ",i.Mv
#    print "  Virial concentration cv = ",i.cv
#    print "     Virial radius Rv/Mpc = ",i.Rv
#    print " "
#    i.otherpars('virial')
#    print "      Scale radius rs/Mpc = ",i.rs
#    print "    Scale rhos/MsunMpc^-3 = ",i.rhos
#    print "     Virial radius Rv/Mpc = ",i.Rv
#    print " "
#    print "Mass difference (%) = ",100*(i.Mv-i.M200)/i.Mv
#    print "rhos difference (%) = ",100*(i.rhos-rhos)/i.rhos
   
   r = 0.29 # Mpc
   rho = i.density(r)
   sigma = i.surface_density(r)
   shear = i.shear(r)
   print "    Density at r=0.29Mpc / MsunMpc^-3 = ",rho
   print "    Surface density at r=0.29Mpc / MsunMpc^-2 = ",sigma
   print "    Shear at r=0.29Mpc / MsunMpc^-2 = ",shear
   
# ======================================================================     
