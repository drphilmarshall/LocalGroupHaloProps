# ======================================================================

import numpy

vb = False

# Age of the Universe:
t0 = 13.73 # Gyr, XXX et al XXXX
t0err = 0.15

# Units etc:
Gyr2sec = 1e9*365.25*24*3600
Mpc2km = 1e6*3.0856e16/1e3
kms2MpcperGyr = Gyr2sec/Mpc2km
# print "To convert km/s to Mpc/Gyr, multiply by",kms2MpcperGyr

G = 4.30e-9 # km^2 s^-2 Mpc Msun^-1

# Solver limits:
tiny = numpy.pi/100.0

# Radial approach:
# Not setting upper limit too close to 2pi was key...
hardxlimits = numpy.array([numpy.pi+tiny,1.9*numpy.pi-tiny])

# Non-radial approach:
# e between 0 and 1 for a closed orbit:
hardelimits = numpy.array([tiny,1.0-tiny])

# ======================================================================
# The timing argument Local Group mass estimate.

def mass(D,vr,vt=None,approach='radial',t0scatter=False):
   
   global fr_constant,ft_constant
   n = len(vr)
   
   # Add scatter to t0?
   if t0scatter:
     t = numpy.random.normal(t0,t0err,n)
   else:
     t = t0

   if (approach == 'radial'):

     fr_constant = (vr * kms2MpcperGyr) * t0 / D
     
     xlimits = numpy.outer(numpy.ones(n),hardxlimits)
     x = solve(fr,xlimits)
     
     # No solution - hide low masses when plotting:
     M = numpy.ones(n)*10.0**9.5
     a = numpy.ones(n)*999.0
     e = numpy.ones(n)
     
     solved = numpy.where(x > hardxlimits[0])
     
     # Sub back for M and a:
     M[solved] = vr[solved]*vr[solved]*D[solved]/(G*(1.0+numpy.cos(x[solved])))
     a[solved] = D[solved]/(1.0 - numpy.cos(x[solved]))

   else:
   # Include tangential speed vt in orbit calculation...

     ft_constant = numpy.sqrt(vr*vr + vt*vt) * kms2MpcperGyr * t0 / D
     
     elimits = numpy.outer(numpy.ones(n),hardelimits)
     e = solve(ft,elimits)
     
     # No solution:
     M = numpy.ones(n)*10.0**10.5
     a = numpy.ones(n)*999.0
     x = numpy.zeros(n)
     esinx = numpy.ones(n)
     ecosx = numpy.ones(n)
     
     solved = numpy.where(e > hardelimits[0])
     
     # Sub back for chi, M and a:
     sinx[solved] = (vr[solved]/vt[solved])*numpy.sqrt(1.0-e[solved]*e[solved])/e[solved]
     ecosx[solved] = e[solved]*numpy.sqrt(1.0-sinx[solved]*sinx[solved])
     x[solved] = numpy.asin(sinx[solved])
     a[solved] = D[solved]/(1.0 - ecosx[solved])
     M[solved] = (vr[solved]*vr[solved]+vt[solved]*vt[solved])*D[solved]/(G*(1.0+ecosx[solved]))

   return M,a,x,e
   
# ----------------------------------------------------------------------
# Function (of chi) to solve in radial approach case.

def fr(x,i):

   global fr_constant
   
   t = fr_constant[i] - numpy.sin(x[i])*(x[i]-numpy.sin(x[i]))/(1.0-numpy.cos(x[i]))**2
   if vb: print "    In fr: x,f = ",x[i],t

   return t
   
# ----------------------------------------------------------------------
# Function (of e) to solve in non-radial approach case.
# UNTESTED

def ft(e,i):

   global ft_constant,vr,vt
   
   sinx = (vr/vt)*numpy.sqrt(1.0-e*e)/e # Closed orbit requires |sinx| < 1
   cosx = numpy.sqrt(1.0-sinx*sinx)
   
   # Warning - no guarantee |sinx| < 1... 
   x = arcsin(sinx) + 2.0*numpy.pi # to ensure [pi:2pi]
   
   t = ft_constant[i] - (x[i] - (e[i]*sinx[i])/(1.0-e[i]*cos(x[i])))*numpy.sqrt((1.0+e[i]*cosx[i])/(1.0-e[i]*cosx[i]))
   if vb: print "    In ft: e,f = ",e[i],t

   return t
   
# ----------------------------------------------------------------------
# Simple bisection rootfinder. Note that x is a vector.

def solve(func,xlimits,tol=1e-5):

   # Initialise arrays:
   x0 = xlimits[:,0]
   x1 = xlimits[:,1]
   n = 1.0*len(x0)
   
   fb = numpy.ones(n)*tol*1000   
   u = numpy.where(fb > tol)   # index marking unsolved systems, here 
                               # used to pick out all systems at start.
   f0 = func(x0,u)
   f1 = func(x1,u)
   df = f1-f0
   xb = numpy.zeros(n)   

   u = numpy.where(f0*f1 > 0.0)  # index marking unsolved systems, here 
                                 # used to flag unsolvable systems.
   fb[u] = 0.0
   xb[u] = x0[u]*f1[u]/df[u] - x1[u]*f0[u]/df[u]
   
   u = numpy.where(numpy.abs(fb) > tol)  # index marking unsolved systems
   if vb: print "  solve: fb = ",fb, "tol = ",tol
   if vb: print "  solve: still working on",u,u[0]
   
   # Now do bisections, updating either x0 or x1 until fb hits tol in 
   # each case - only have to calculate for the u-indexed systems:   
   i = 0
   while len(u[0]) > 0:
 
      i += 1
      if vb: print "  solve: iteration",i,", unsolved count = ",len(u[0])
 
      # Evaluate function at current limits, and reset unfinished index:
      f0[u] = func(x0,u)
      f1[u] = func(x1,u)
      df[u] = f1[u] - f0[u]
      xb[u] = x0[u]*f1[u]/df[u] - x1[u]*f0[u]/df[u]
      fb[u] = func(xb,u)
      
      un = (numpy.abs(fb) > tol)
      u = numpy.where(un)
      if vb: print "   solve: xb = ",xb
      if vb: print "   solve: fb = ",fb, "tol = ",tol
      if vb: print "   solve: still working on",u,u[0]
      
      # Update limits to better bracket the roots:
      m0 = numpy.where((un) & (numpy.sign(fb) == numpy.sign(f0)))
      m1 = numpy.where((un) & (numpy.sign(fb) == numpy.sign(f1)))
      x0[m0] = xb[m0]
      x1[m1] = xb[m1]
      if vb: print "   solve: updated x0,x1 = ",x0,x1

   return xb


# ======================================================================     

# Testing:

if __name__ == '__main__':

   vr = numpy.array([-120.0,-130.0,-140.0])
   D = numpy.array([700.0,800.0,900.0])/1000.0
   vb = False
   
   print "Radial velocity vr/kms^-1 = ",vr
   print "Distance D/Mpc = ",D
   
   M,a,chi = mass(vr,D,t0scatter=True)
   
   print "Mass estimates M/Msun = ",M
   print "Semi-major axes a/Mpc = ",a
   print "Eccentric anomalies chi/pi = ",chi/numpy.pi
