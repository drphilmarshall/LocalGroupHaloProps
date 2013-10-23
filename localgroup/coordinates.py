# ======================================================================

from numpy import sin,cos,arcsin,arccos,pi,abs,sqrt

deg2rad = pi/180.0

# North Galactic pole:
alphaG = 192.85948 * deg2rad
deltaG =  27.12825 * deg2rad
# Galactic longitude of the ascending node of the Galactic plane:
lOmega = 32.93192 * deg2rad

# Numerical tolerance:
tol = 1e-8

# ======================================================================
# Rotate sky positions from equatorial to galactic coordinates.
# From Poleski 2013.

def equatorial_to_galactic(RA,DEC):
    
    alpha = RA*deg2rad
    delta = DEC*deg2rad
    
    b = arcsin(cos(delta)*cos(deltaG)*cos(alpha - alphaG) + sin(delta)*sin(deltaG))
    l = lOmega - arcsin(cos(delta)*sin(alpha - alphaG)/cos(b))
    
    assert abs(cos(lOmega - l)*cos(b) - sin(delta)*cos(deltaG)-cos(delta)*sin(deltaG)*cos(alpha-alphaG)) < tol
    
    return l/deg2rad,b/deg2rad
    
# ----------------------------------------------------------------------
# Rotate proper motions from equatorial to galactic coordinates.
# From Poleski 2013.

def equatorial_to_galactic_proper_motion(vW,vN,RA,DEC):
    
    alpha = RA*deg2rad
    delta = DEC*deg2rad
    
    C1 = sin(deltaG)*cos(delta) - cos(deltaG)*sin(delta)*cos(alpha-alphaG)
    C2 = cos(deltaG)*sin(alpha-alphaG)
    cosb = sqrt(C1*C1+C2*C2)
    
# When working with v's, we can drop the extra cosines, I think...
#     vl = (C1*vW*cos(delta) + C2*vN           )/(cosb*cosb)
#     vb = (C1*vN            - C2*vW*cos(delta))/cosb

    vl = (C1*vW + C2*vN)/cosb
    vb = (C1*vN - C2*vW)/cosb

    return vl,vb
    
# ----------------------------------------------------------------------
# Generic transformation between spherical and cartesian coordinates.
# For galactic coords, use l,b instead of RA,DEC, and so on.

def spherical_to_cartesian(RA,DEC,D,v_west,v_north,v_r):

    # Heliocentric Cartesian positions (Mpc):
    delta = DEC*deg2rad
    alpha = RA*deg2rad
    x = D*cos(delta)*cos(alpha)
    y = D*cos(delta)*sin(alpha)
    z = D*sin(delta)
     
    # Heliocentric Cartesian velocities (km/s):
    vx = (x/D)*v_r - (z/D)*v_north*cos(alpha) - v_west*sin(alpha) 
    vy = (y/D)*v_r - (z/D)*v_north*sin(alpha) + v_west*cos(alpha) 
    vz = (z/D)*v_r + cos(delta)*v_north

    return x,y,z,vx,vy,vz

# ======================================================================
