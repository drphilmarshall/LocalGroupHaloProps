# ======================================================================

from numpy import mean,sin,cos,arcsin,arccos,arctan2,pi,abs,sqrt,sign,atleast_1d,arctan,std, array
# from pytpm import tpm, convert
import numpy as np
import math

deg2rad = pi/180.0
arcsec2rad = deg2rad/3600.0
muaspyrMpc2kmps = 4.7404

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
    b = atleast_1d(b)

    sindl = cos(delta)*sin(alpha-alphaG)/cos(b)
    cosdl = (sin(delta)*cos(deltaG) - cos(delta)*sin(deltaG)*cos(alpha-alphaG))/cos(b)
    dl = arctan2(cosdl,sindl)
    l = atleast_1d(lOmega + dl)
    l[l < 0.0] += 2.0*pi


    return l/deg2rad,b/deg2rad

# ----------------------------------------------------------------------
# Rotate proper motions from equatorial to galactic coordinates.
# From Poleski 2013.

def equatorial_to_galactic_proper_motion(mu_w, mu_n,RA,DEC):

    mu_alpha = mu_w
    mu_delta = mu_n
    alpha = RA*deg2rad
    delta = DEC*deg2rad

    C1 = sin(deltaG)*cos(delta) - cos(deltaG)*sin(delta)*cos(alpha-alphaG)
    C2 = cos(deltaG)*sin(alpha-alphaG)
    cosb = sqrt(C1*C1+C2*C2)

    mu_l = (C1*mu_alpha + C2*mu_delta)/cosb
    mu_b = (C1*mu_delta - C2*mu_alpha)/cosb


    return mu_l,mu_b

# ----------------------------------------------------------------------
# Generic transformation between spherical and cartesian coordinates.
# For galactic coords, use l,b instead of RA,DEC, and so on.
# UNITS:  RA, DEC should be in degrees.  D is in Mpc.
# mu_[west|north] in mas/yr, v_r in km/s
# x,y,z in Mpc.  vx, vy, vz, in km/s

def spherical_to_cartesian(RA,DEC,D,mu_west,mu_north,v_r,halo):

    # Heliocentric Cartesian positions (Mpc):
    delta = DEC*deg2rad
    alpha = RA*deg2rad
    x = D*cos(delta)*cos(alpha)
    y = D*cos(delta)*sin(alpha)
    z = D*sin(delta)

    # Convert mu to v
    if not halo=='M31':
        v_west = -mu_west*muaspyrMpc2kmps*D
        v_north = mu_north*muaspyrMpc2kmps*D
    else:
        v_west = -mu_west
        v_north = mu_north
    print 'v_west:  %f +/- %f'%(mean(v_west),std(v_west)) 
    print 'v_north: %f +/- %f'%(mean(v_north),std(v_north))

    return v_west, v_north

# ----------------------------------------------------------------------
# Transformation from Heliocentric Galactic Cartesian system to Galactocentric Cartesian.
# NOTE: The different definitions of the axis directions.

def heliocentric_galactic_cartesian_to_galactocentric_cartesian(xh,yh,zh,vxh,vyh,vzh,R0, VX, V0, VZ):

    xg = xh - R0
    yg = yh
    zg = zh
    vxg = vxh - VX
    vyg = vyh - V0
    vzg = vzh - VZ

    return xg,yg,zg,vxg,vyg,vzg

def unit_vec(vec):
    vec = np.asanyarray(vec)
    vec /= math.sqrt((vec*vec).sum())
    return vec

def get_swn_basis(l, b):
    l *= math.pi/180.0
    b *= math.pi/180.0
    sys = unit_vec([math.cos(l)*math.cos(b), math.sin(l)*math.cos(b), math.sin(b)]) #los
    west = unit_vec(np.cross(sys, np.array([0,0.0,1.0])))
    north = unit_vec(np.cross(west, sys))
    return np.vstack([sys, west, north])

def get_cartesians(ra, dec, d, vw, vn, vr, R0, VX, VY, VZ):
    gal_vel_swn = np.array([vr, -vw, vn])
    sun_vel_gc = np.array([VX, VY, VZ])
    gal_swn_basis = get_swn_basis(ra, dec)
    vx,vy,vz = np.dot(gal_vel_swn, gal_swn_basis) - sun_vel_gc
    x, y, z = gal_swn_basis[0]*d - np.array([R0, 0, 0])
    return x, y, z, vx, vy, vz

# ----------------------------------------------------------------------
# Transforms from Heliocentric equatorial spherical coordinates (ra, dec, etc...) to
# galactocentric cartesian coordinates (x,y,z, etc).
# NOTE:  Add optional argument to take care of internal rotation.
#def heliocentric_equatorial_spherical_to_galactocentric_cartesian(ra, dec, d, mu_w, mu_n, v_r, dvrot_w, dvrot_n, R0=0.00829, VX=-11, V0=-239, VZ=-7, halo='M31'):

def heliocentric_equatorial_spherical_to_galactocentric_cartesian(ra, dec, d, mu_w, mu_n, v_r, R0, VX, V0, VZ, halo='M31'):
    print 'Inside method: heliocentric_equatorial_spherical_to_galactocentric_cartesian'
    l,b = equatorial_to_galactic(ra,dec)
    print 'l = ',l
    print 'b = ',b

    if halo=='M31':
        mu_l = mu_w
        mu_b = mu_n
        mu_l,mu_b = equatorial_to_galactic_proper_motion(-mu_w,mu_n,ra,dec)
        mu_l = -mu_l
    elif halo=='M33':
        mu_l = mu_w
        mu_b = mu_n
        mu_l,mu_b = equatorial_to_galactic_proper_motion(-mu_w,mu_n,ra,dec)
        mu_l = -mu_l
    elif halo=='LMC':
        mu_l,mu_b = equatorial_to_galactic_proper_motion(-mu_w,mu_n,ra,dec)
        mu_l = -mu_l
    else:
        raise ValueError('unknown galaxy name')


    #mu_l,mu_b = equatorial_to_galactic_proper_motion(-mu_w,mu_n,ra,dec)
    print 'mu_l = %f +/- %f'%(mean(mu_l), std(mu_l))
    print 'mu_b = %f +/- %f'%(mean(mu_b), std(mu_b))
    v_west, v_north = spherical_to_cartesian(l, b, d, mu_l, mu_b, v_r, halo)
    x = []
    y = []
    z = []
    vx = []
    vy = []
    vz = []
    for i in range(len(ra)):
        xx, yy, zz, vxx, vyy, vzz = get_cartesians(l[i], b[i], d[i], v_west[i], v_north[i], v_r[i], R0[i], VX[i], V0[i], VZ[i])
        x.append(xx)
        y.append(yy)
        z.append(zz)
        vx.append(vxx)
        vy.append(vyy)
        vz.append(vzz)


    return np.array(x),np.array(y),np.array(z),np.array(vx),np.array(vy),np.array(vz)

    return x,y,z,vx,vy,vz

# ======================================================================

# Test on Hipparcos stars, from XHIP catalog
# http://www.astrostudio.org/xhipreadme.html
# accessible by Vizier at http://vizier.u-strasbg.fr/viz-bin/VizieR-4

# # HIP   RAJ2000       DEJ2000         pmRA    pmDE    GLon            GLat          Dist    pmGLon   pmGLat    X      Y      Z      RV     U       V      W
# #       deg           deg             mas/yr  mas/yr  deg             deg           pc      mas/yr   mas/yr    pc     pc     pc     km/s   km/s    km/s   km/s
# 950     2.93286720    -35.13339602    169.76  114.63  347.15730758    -78.33861529  21.28   -23.51  -203.48     4.2   -1.0   -20.8  -21.0  -24.3     3.1   16.5
# 1013    3.16145707    -26.85334224    266.51  115.84   32.54473512    -81.35606282  70.48   134.25  -257.73     8.9    5.7   -69.7    5.2  -95.2    -7.6  -18.1
# 1092    3.40754870     80.66533914    251.57  182.95  121.32780663     17.91155946  19.56   276.69   142.14    -9.7   15.9     6.0  -13.4  -13.2   -27.7    8.4
# 1297    4.06050131     35.17799681    156.98  131.63  114.85499911    -27.14115851  66.49   175.27   106.07   -24.9   53.7   -30.3  -18.1  -49.8   -24.0   38.0
# 1349    4.22329979    -52.65159218    312.90  180.90  314.78571961    -63.67436085  22.59  -243.84  -266.78     7.1   -7.1   -20.3  -11.7  -40.2     3.5   -2.2
# 1386    4.33276003     29.18196398    689.56  427.58  114.04242624    -33.10037842  23.12   748.30   313.62    -7.9   17.7   -12.6  -10.4  -79.0   -24.2   34.5
# 1475    4.58559072     44.02195596   2887.50  408.90  116.66943324    -18.44699193   3.59  2916.27    15.31    -1.5    3.0    -1.1   11.6  -49.3   -12.3   -3.4

# ======================================================================

if __name__ == '__main__':

    #RA,DEC = 2.93286720,-35.13339602
    RA, DEC = .14287, -53.0976
    #l,b = 347.15730758,-78.33861529
    l,b = 96.13532, -60.538212
    ll,bb = equatorial_to_galactic(RA,DEC)
    print "l,b = ",ll,bb
    print "  cf XHIP values: ",l,b

    #mu_a,mu_d = 169.76,114.63
    mu_a, mu_d = 62.20, -15.30
    mu_l,mu_b = equatorial_to_galactic_proper_motion(mu_a,mu_d,RA,DEC)
    print "mu_l,mu_b = ",mu_l,mu_b
    #mu_ll, mu_bb = -23.51, -203.48
    mu_ll, mu_bb = 50.83, -38.97
    print "  cf XHIP values: ", mu_ll, mu_bb

    #D = 21.28
    D = 201.17
    #v_r = -21.0
    v_r = 21.61
    d_vn = 0
    d_vw = 0
    K = 4.7404e-3
    v_l = K*mu_l*D
    v_b = K*mu_b*D

    x,y,z,vx,vy,vz = spherical_to_cartesian(l,b,D,mu_l/1000,mu_b/1000,v_r,d_vn, d_vw)

    print "x,y,z = ",x,y,z
    print "  cf XHIP values: 4.2, -1.0, -20.8"
    print "vx,vy,vz = ",vx,vy,vz
    print "  cf XHIP values: -24.3, 3.1, 16.5"

    print "OK FINE"

    print "Now check galactic center (l,b) = 0,0:"
    xh,yh,zh,vxh,vyh,vzh = spherical_to_cartesian(0,0,0.0085, -5460.,0.,0., d_vn, d_vw)
    print "Heliocentric Galactic Cartesian location of G.C.: ", xh,yh,zh,vxh,vyh,vzh

    x,y,z,vx,vy,vz = heliocentric_galactic_cartesian_to_galactocentric_cartesian(xh, yh, zh, vxh, vyh, vzh)
    print "Galactocentric cartesian coordinates of G.C. (should be all zero): ",x,y,z,vx,vy,vz

    x,y,z,vx,vy,vz = heliocentric_equatorial_spherical_to_galactocentric_cartesian((17.0/24.+42./(24.0*60.0)+29.319/(24.*3600.))*360., -28-59./60-18.54/3600., 0.0085, 2.7, -5.6, 0.0,0,0)
    print "Galactocentric cartesian coordinates of G.C. (should be all zero): ",x,y,z,vx,vy,vz

    print ""

# l,b =  347.157307577 -78.3386152875
#   cf XHIP values:  347.15730758 -78.33861529
# mu_l,mu_b =  -23.5134553445 -203.483689562
#   cf XHIP values: -23.51,-203.48
# x,y,z =  4.19366817246 -0.95606468303 -20.8407650431
#   cf XHIP values: 4.2, -1.0, -20.8
# vx,vy,vz =  -24.2656956924 3.09925289521 16.4175661185
#   cf XHIP values: -24.3, 3.1, 16.5
# ======================================================================
