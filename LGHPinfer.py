#!/usr/bin/env python
# ======================================================================

# Globally useful modules, imported here and then accessible by all
# functions in this file:

import numpy,atpy,sys,getopt,string,subprocess,pyfits

import localgroup

# It would be nice to read in tables from FITS, but I cannot write them 
# out in FITS format yet :-/  Stick with ascii now:
# tabletype = 'fits'
tabletype = 'ascii'

rad2deg = 180.0/numpy.pi

# ======================================================================

def LGHPinfer(argv):

  """
  NAME
    LGHPinfer.py

  PURPOSE
    Read in catalogs of paired halos in known format, and join to make a  list
    of local group "triplets". Compute some new quantities including the
    importance of each sample given various observational constraints. Output
    a plain text catalog for use by CornerPlotter.

  COMMENTS
    If no input files are provided, 3 default files will be read in, assumed to 
    contain triplet halo samples from the prior. These files are called 
    LG_halos_M[W,31,33].fits - although they may not be in FITS format...

  USAGE
    LGHPinfer.py [flags] [options] pairs.dat

  FLAGS
          -h               Print this message [0]
          -v               Be verbose

  INPUTS
          MW_halos_01.txt     Busha format plain text catalog
          M31_halos_01.txt    Busha format plain text catalog
          M33_halos_01.txt    Busha format plain text catalog

  OPTIONAL INPUTS
          
  OUTPUTS
          stdout       Useful information
          *.cpt        Catalogs for plotters to read

  EXAMPLES
  
     LGHPinfer.py -c constraints.txt  *halos*.txt

  BUGS


  HISTORY
    2011-10-21 started Marshall & Busha (Oxford)
    2012-03-15 added kinematic constraints Marshall (Oxford)
    2013-10-01 githubbed Marshall & Busha (KIPAC)
  """

  # --------------------------------------------------------------------

  try:
      opts, args = getopt.getopt(argv[1:], "hvr:", ["help","verbose","vt=","vr=","M1="])
  except getopt.GetoptError, err:
      # print help information and exit:
      print str(err) # will print something like "option -a not recognized"
      print LGHPinfer.__doc__
      sys.exit(2)

  vb = False
  headstart = False
  constrainMMW = False
  confile = None  
  for o,a in opts:
    if o == "-v":
        vb = True
    elif o in ("-h", "--help"):
        print LGHPinfer.__doc__
        return
    elif o in ("-MW"):
        bits = a.split(',')
        if len(bits) != 2:
          print "ERROR: supply MW constraint as 'M1bar,M1err'"
          return
        MMWbar = bits[0]
        MMWerr = bits[1]
        constrainMMW = True  
    elif o in ("-c", "--constraints"):
        confile = a
    else:
        assert False, "unhandled option"

  if (len(args) > 0) and (len(args) % 3 == 0):
    inputs = args
  elif (len(args) == 0):
    headstart = True
  else:
    print LGHPinfer.__doc__
    return
  
  # If FITS tables of halo triplets exist, read them in (fast):
  
  if headstart:
  
     if vb: print "Working from existing FITS tables of halos..."
  
  # --------------------------------------------------------------------
  # PART 1: reading in raw ascii data, and selecting triplets. Output 
  # of this part is 3 FITS files, one for each galaxy, with all relevant 
  # quantities.
  
  else:
     # Parse list of input files, and read in data as appended tables:

     started = {'MW':0,'M31':0,'M33':0}

     for input in inputs:

       pieces = string.split(string.split(input,'.')[0],'_')
       object = pieces[0]
       box = pieces[2]
       epoch = pieces[3]

       # Read in table and add a column of strings for the consuelo run no:
       if vb: print "Reading data from "+input
       t = atpy.Table(input,type='ascii')
       consuelo_box = numpy.zeros(len(t),dtype='|S4')
       consuelo_box.fill(box)
       t.add_column('run',consuelo_box)
       consuelo_epoch = numpy.zeros(len(t),dtype='|S6')
       consuelo_epoch.fill(epoch)
       t.add_column('epoch',consuelo_epoch)

       if not started[object]:
       # Rename table to new one!
          if (object == 'MW'):
            MW = t
          elif (object == 'M31'):
            M31 = t
          elif (object == 'M33'):
            M33 = t
          else:
            print "Error: Unrecognised object: "+object
            return
          started[object] = 1

       else:    
       # Append table to one that has begun:
          if (object == 'MW'):
            MW.append(t)
          elif (object == 'M31'):
            M31.append(t)
          elif (object == 'M33'):
            M33.append(t)
          else:
            print "Error: Unrecognised object: "+object
            return

     # Now have all tables read in and combined into 3, one of r each object.

     if vb:
       print "Read in",len(MW),"MW halos, ",len(M31),"M31 halos, and",len(M33),"M33 halos"
     if (len(MW) != len(M31)):
       print "Error: mismatched dataset sizes: ",len(MW),"MW halos cf",len(M31),"M31 halos"
       return

     # Tables are matched via IDs, so sort by M33 ID, run, and MW_ID and reverse
     # to get high idM33 values at the top. Reversing idM33 vals 
     # needs to be done first!

     sort_to_align(MW)

     sort_to_align(M31)

     # Check that arrays line up:
     fail = numpy.sum(MW.idMW - M31.idMW)
     if fail:
       print "Error: IDs in data files not consistent"
       print "  MW.idMW = ",MW.idMW
       print "  M31.idMW = ",M31.idMW
       return

     # Now sort M33 to match:

     sort_to_align(M33)

     # Can now make index array for M33 triplets with numpy.where, and use it 
     # on MW and M31 tables: this is basically the first N samples in all 3
     # tables. If M33 is not being used, can use all M31 and MW halos anyway.

     index = numpy.where(MW.idM33 >= -2)
     N_pairs = len(index[0])
     indexM33 = numpy.where(MW.idM33 >= 0)
     N_triplets = len(indexM33[0])
     if vb: print "Found",N_pairs,"MW-M31 pairs, and",N_triplets,"LG triplets"

     # Check that arrays line up:
     if N_triplets != len(M33):
       print "Warning: numbers of M33s in datafiles are inconsistent:"
       print "  len(MW.idM33[indexM33]), len(M33.idM33) = ",len(MW.idM33[indexM33]), len(M3.idM33)
       print "  MW.idM33[indexM33] = ",MW.idM33[indexM33]
       print "  M33.idM33 = ",M33.idM33
       return

     fail = numpy.sum(MW.idM33[indexM33] - M33.idM33)
     if fail:
       print "Warning: IDs in data files are inconsistent:"
       print "  MW.idM33[indexM33] = ",MW.idM33[indexM33]
       print "  M33.idM33 = ",M33.idM33
       return

     #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     # Shift all halos to galactocentric coordinates, subtracting off MW
     # position and velocity - and rescale to include Hubble flow. Note 
     # that distances in catalog are in Mpc!:

     H0 = 70.0
     phase_space_shift_and_flow(MW,M31,M33,indexM33,H0)

     # Add columns with galactocentric distances D and galactocentric 
     # radial velocity v (possibly slow):

     distances_and_velocities(M31,M33) 

     # Add columns with M200 and c200 to complement Mvir and cvir

     massconvert(MW)
     massconvert(M31)
     massconvert(M33)

     #  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

     # Make a unit importance array, and write out cpt files with just the 
     # prior samples in them. Logs of mass are taken just before writing.

     w = numpy.ones(len(MW))
     w_M33 = numpy.ones(len(M33))

     # BUG: Some of Michael's halos are hundreds of Mpc away - some sort of 
     # wraparound error. Re-impose the isolation selection via the weights:

     wrapped = numpy.where(M31.D > 1.5)
     w[wrapped] = 0.0
     wrapped = numpy.where(M33.D > 1.5)
     w_M33[wrapped] = 0.0

     # BUG: Occasionally a triplet has M31 = M33. Remove these!
     identical = numpy.where(M31.idM31 == M31.idM33)
     w[identical] = 0.0
     w_M33[identical] = 0.0

     # At this point, write out plain text catalogs for visualising prior - 
     # random n samples:

     # n = 2000
     # 
     # filename = 'fig1_pairs_prior.cpt'
     # write_cpt_file(index,MW,M31,M33,w,n,filename)
     # if vb: print "Written",n,"lines to",filename
     # 
     # filename = 'fig1_triplets_prior.cpt'
     # write_cpt_file(indexM33,MW,M31,M33,w_M33,n,filename)
     # if vb: print "Written",n,"lines to",filename
     
     
     # Write out pairs and triplets for future use (save having to read 
     # in an enormous file 
     # again). Weights are added at this point:  

     write_fits_files('pairs',index,MW,M31,M33,w,vb)
     write_fits_files('triplets',indexM33,MW,M31,M33,w_M33,vb)
   
     # OK, now have 3+2 tables, each with a weight column, 
     # identical for now.

     # # What's the probability that M_MW is > M_M31?
     # P_pairs = probgt(index,MW,M31,w)
     # P_triplets = probgt(indexM33,MW,M31,w)
     # print "Probability that M_MW > M_M31 is",P_pairs*100.0,"% for pairs,"
     # print "  and",P_triplets*100.0,"% for triplets"

  # End of Part 1.
  
  # --------------------------------------------------------------------

  # Part 2: compute likelihood of observational data, and apply to each 
  # triplet as an importance for drawing inferences. First need some 
  # additional quantities, beyond D and vr.
  
  
  # Read in data tables afresh, from pre-prepared files.
  
  input = 'LG_triplets_MW.fits'
  MWt = atpy.Table(input,type=tabletype,name='MW')
  input = 'LG_triplets_M31.fits'
  M31t = atpy.Table(input,type=tabletype,name='M31')
  input = 'LG_triplets_M33.fits'
  M33t = atpy.Table(input,type=tabletype,name='M33')
  if vb: print "Read in",len(MWt),"MW halos,",len(M31t),"M31 halos, and",len(M33t),"M33 halos in triplets"

#   input = 'LG_pairs_MW.fits'
#   MWp = atpy.Table(input,type=tabletype,name='MW')
#   input = 'LG_pairs_M31.fits'
#   M31p = atpy.Table(input,type=tabletype,name='M31')
#   if vb: print "Read in",len(MWp),"MW halos and",len(M31p),"M31 halos in pairs"

  # Just make sure they all agree:
  fail = ((2*len(MWt)-len(M31t)-len(M33t)) != 0)
  if fail:
    print "Warning: triplet FITS tables have inconsistent lengths:"
    print "  len(MWt) = ",len(MWt)
    print "  len(M31t) = ",len(M31t)
    print "  len(M33t) = ",len(M33t)
    return
#   fail = ((len(MWp)-len(M31p)) != 0)
#   if fail:
#     print "Warning: triplet FITS tables have inconsistent lengths:"
#     print "  len(MWp) = ",len(MWp)
#     print "  len(M31p) = ",len(M31p)
#     return
#   
#   pindex = numpy.where(MWp.idM31 >= 0)
#   np = len(pindex[0])
  tindex = numpy.where(MWt.idM33 >= 0)
  nt = len(tindex[0])
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  # Add column with angsep to M33: 

  angular_separation_on_sky(M31t,M33t) 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  # Add column with timing argument mass, to M31 tables: 

#   timing_argument(M31p) 
  timing_argument(M31t) 

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  
  # Start out with prior weight (all ones except a few zeros):
#   wp0 = MWp.w * M31p.w
  wt0 = MWt.w * M31t.w * M33t.w

  # Read in constraints as an instance of the constraints class, which 
  # contains the relevant observational data:
    
  obs = LGconstraints.bundle()
  
  # Calculate importances based on input MW and M31 constraints.
  
  # Triplets:

#   # Compute MW mass constraint, and save the posterior:
#   wt1 = wt0*obs.likelihood(MWt,'Mvir')
#   filename = 'LG_triplets_MW-Mvir.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt1,nt,filename,vb)
  
  # Compute M31 distance constraint, and save the posterior:
  wt2 = wt0*obs.likelihood(M31t,'D')
#   filename = 'LG_triplets_M31_D.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt2,nt,filename,vb)
  
  # Compute M31 radial velocity constraint, and save the posterior:
  wt3 = wt0*obs.likelihood(M31t,'vr')
#   filename = 'LG_triplets_M31_vr.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt3,nt,filename,vb)
  
  # Compute M33 distance constraint, and save the posterior:
  wt4 = wt0*obs.likelihood(M33t,'D')
#   filename = 'LG_triplets_M33_D.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt4,nt,filename,vb)
  
  # Compute M33 radial velocity constraint, and save the posterior:
  wt5 = wt0*obs.likelihood(M33t,'vr')
#   filename = 'LG_triplets_M33_vr.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt5,nt,filename,vb)
  
  # Combine all kinematic constraints, and save the posterior:
  wt6 = wt2*wt3*wt4*wt5
  filename = 'LG_triplets_M33+M31_D+vr.cpt'
  write_cpt_file(tindex,MWt,M31t,M33t,wt6,nt,filename,vb)
  
#   # Combine just M31 kinematic constraints, and save the posterior:
#   wt7 = wt2*wt3
#   filename = 'LG_triplets_M31_D+vr.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt7,nt,filename,vb)
#   
#   # Combine all constraints, and save the posterior:
#   wt8 = wt1*wt2*wt3*wt4*wt5
#   filename = 'LG_triplets_M33+M31_D+vr_MW-Mvir.cpt'
#   write_cpt_file(tindex,MWt,M31t,M33t,wt8,nt,filename,vb)
  
  # Write out prior for completeness:
  filename = 'LG_triplets_prior.cpt'
  write_cpt_file(tindex,MWt,M31t,M33t,wt0,nt,filename,vb)
  
#   # Pairs:
# 
#   # Compute MW mass constraint, and save the posterior:
#   wp1 = wp0*obs.likelihood(MWp,'Mvir')
#   filename = 'LG_pairs_MW_Mvir.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp1,np,filename,vb)
#   
#   # Compute M31 distance constraint, and save the posterior:
#   wp2 = wp0*obs.likelihood(M31p,'D')
#   filename = 'LG_pairs_M31_D.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp2,np,filename,vb)
#   
#   # Compute M31 radial velocity constraint, and save the posterior:
#   wp3 = wp0*obs.likelihood(M31p,'vr')
#   filename = 'LG_pairs_M31_vr.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp3,np,filename,vb)
#   
#   # Combine both kinematic constraints, and save the posterior:
#   wp6 = wp2*wp3
#   filename = 'LG_pairs_M31_D+vr.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp6,np,filename,vb)
#   
#   # Combine all constraints, and save the posterior:
#   wp8 = wp1*wp2*wp3
#   filename = 'LG_pairs_M31_D+vr_MW-Mvir.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp8,np,filename,vb)
#   
#   # Write out prior for completeness:
#   filename = 'LG_pairs_prior.cpt'
#   write_cpt_file(pindex,MWp,M31p,M31p,wp0,np,filename,vb)
#   
#   
#   # Write out FITS files for Risa - all constraints except MW mass.
#   
#   write_fits_files('pairs_constrained-by-M31-D+vr',pindex,MWp,M31p,M31p,wp6,vb)
#   write_fits_files('triplets_constrained-by-M33+M31_D+vr',tindex,MWt,M31t,M33t,wt6,vb)
#   
  # --------------------------------------------------------------------

  return

# ======================================================================

def sort_to_align(t):

  t.sort(['idM33','run','idMW'])
  #   print "Top 5 IDs: ",t.idM33[0:5]
  reverse = numpy.arange(len(t.idM33),0,-1)
  t.add_column('reverse',reverse)
  t.sort(['reverse'])
  t.remove_columns(['reverse'])
  #   print "After reverse sorting, Top 5 IDs: ",t.idM33[0:5]
  
  return
  
# ======================================================================
  
def phase_space_shift_and_flow(MW,M31,M33,i,H0):

  M31.x  -= MW.x
  M31.y  -= MW.y
  M31.z  -= MW.z
  M33.x  -= MW.x[i]
  M33.y  -= MW.y[i]
  M33.z  -= MW.z[i]

  MW.x   -= MW.x
  MW.y   -= MW.y
  MW.z   -= MW.z
  
  M31.vx -= MW.vx
  M31.vy -= MW.vy
  M31.vz -= MW.vz
  M33.vx -= MW.vx[i]
  M33.vy -= MW.vy[i]
  M33.vz -= MW.vz[i]

  MW.vx  -= MW.vx
  MW.vy  -= MW.vy
  MW.vz  -= MW.vz
  
  M31.vx += H0*M31.x
  M31.vy += H0*M31.y
  M31.vz += H0*M31.z
  M33.vx += H0*M33.x
  M33.vy += H0*M33.y
  M33.vz += H0*M33.z
  
  return
  
# ======================================================================
  
def distances_and_velocities(M31,M33):
  
  D = numpy.sqrt(M31.x*M31.x + M31.y*M31.y + M31.z*M31.z)
  ii = numpy.where(D <= 0.0)
  if len(ii[0]) > 0:
    print "Error: M31 distance is zero or negative:"
    for i in ii[0]:
      print M31[i]
#       D[i] = 0.75
    sys.exit()
  M31.add_column('D',D)
  
  D = numpy.sqrt(M33.x*M33.x + M33.y*M33.y + M33.z*M33.z)
  ii = numpy.where(D <= 0.0)
  if len(ii[0]) > 0:
    print "Error: M33 distance is zero or negative:"
    for i in ii[0]:
      print M33[i]
#       D[i] = 0.75
    sys.exit()
  M33.add_column('D',D)

  # Project v onto radial vector to get vr, and onto the sphere to 
  # get vt:
  v = numpy.sqrt(M31.vx*M31.vx + M31.vy*M31.vy + M31.vz*M31.vz)
  vr = (M31.vx*M31.x + M31.vy*M31.y + M31.vz*M31.z)/M31.D
  costheta = vr/v
  sintheta = numpy.sqrt(1.0-costheta*costheta)
  vt = v*sintheta
  M31.add_column('vr',vr)
  M31.add_column('vt',vt)
  
  v = numpy.sqrt(M33.vx*M33.vx + M33.vy*M33.vy + M33.vz*M33.vz)
  vr = (M33.vx*M33.x + M33.vy*M33.y + M33.vz*M33.z)/M33.D
  costheta = vr/v
  sintheta = numpy.sqrt(1.0-costheta*costheta)
  vt = v*sintheta
  M33.add_column('vr',vr)
  M33.add_column('vt',vt)
  
  return
  
# ======================================================================
# Compute angular separation on sky - ignore Earth's position 
# relative to MW centre for now... This will need including though!
# Just dot the position vectors together for now.

def angular_separation_on_sky(M31,M33):
    
  cos_angsep = (M31.x*M33.x + M31.y*M33.y + M31.z*M33.z)/(M31.D*M33.D)
  angsep = numpy.arccos(cos_angsep)*rad2deg

  # Note - pay attention to weights to avoid identical twins:
  index = numpy.where(M31.w == 0.0)
  angsep[index] = 0.0
  
  M33.add_column('angsep',angsep)
  
  return
  
# ======================================================================
# Convert virial masses to M200c, to better compare with Li & White.

def massconvert(table):
    
   NFW.vb = True
   
   # If we had z and c for each halo:
   # h = NFW.halo(table.Mvir,table.z,c=table.cvir,kind='virial')   
   
   # Placeholder:
   z = 0.0
   h = NFW.halo(table.Mvir,z,kind='virial')   
   
   h.massconvert('virial')
  
   table.add_column('M200',h.M200)
   table.add_column('c200',h.c200)
  
   return
  
# ======================================================================
# Compute timing argument LG mass, M_TA

def timing_argument(M31):
      
  TimingArgument.vb = False
  M,a,chi,e = TimingArgument.mass(M31.D,M31.vr,approach='radial',t0scatter=True)
  
  M31.add_column('M_TA',M)
  
  return
  
# ======================================================================
# Return the probability that A is greater than B:

def probgt(index,A,B,w):

  # Check that weights sum to 1.0:
  ww = w[:]
  norm = numpy.sum(w[index])
  ww[index] = w[index] / norm
  
  # Now sum probabilities:
  return numpy.sum(ww[numpy.where(A[index] > B[index])])
  
# ======================================================================

def write_cpt_file(index,MW,M31,M33,w,n,filename,vb):

  n_given = len(index[0])
  n_required = n
  if n_required > n_given:  n_required = n_given

# First define labels and ranges:

  logM_MW_label           = '$\log_{10} M_{\\rm MW} / M_{\odot}$,  '
  logM_MW_range           = '10.0,14.0,    '

  logM_M31_label          = '$\log_{10} M_{\\rm M31} / M_{\odot}$,  '
  logM_M31_range          = '10.0,14.0,    '

  logM_M33_label          = '$\log_{10} M_{\\rm M33} / M_{\odot}$,  '
  logM_M33_range          = '10.0,14.0,    '

  logM_LG_pair_label      = '$\log_{10} M\prime_{\\rm LG} / M_{\odot}$,  '
  logM_LG_pair_range      = '10.0,14.0,    '

  logM_LG_triplet_label   = '$\log_{10} M_{\\rm LG} / M_{\odot}$,  '
  logM_LG_triplet_range   = '10.0,14.0,    '

  logMratio_label         = '$\log_{10} M_{\\rm M31} / M_{\\rm MW}$,  '
  logMratio_range         = '-3.0,3.0,     '

  D_M31_label             = '$D_{\\rm M31} / {\\rm kpc}$,  '
  D_M31_range             = '550,950,      '

  D_M33_label             = '$D_{\\rm M33} / {\\rm kpc}$,  '
  D_M33_range             = '550,950,      '

  vr_M31_label            = '$v^{\\rm rad}_{\\rm M31} / {\\rm km s}^{-1}$,  '
  vr_M31_range            = '-400,400,     '

  vt_M31_label            = '$v^{\\rm tan}_{\\rm M31} / {\\rm km s}^{-1}$,  '
  vt_M31_range            = '0,400,     '

  vr_M33_label            = '$v^{\\rm rad}_{\\rm M33} / {\\rm km s}^{-1}$,  '
  vr_M33_range            = '-400,400,     '
  
  vt_M33_label            = '$v^{\\rm tan}_{\\rm M33} / {\\rm km s}^{-1}$,  '
  vt_M33_range            = '0,400,     '

  angsep_label            = '$\\Delta\\theta / {\\rm deg}$,  '
  angsep_range            = '0,60,     '
  
  logM_TA_label           = '$\log_{10} M_{\\rm TA} / M_{\odot}$,  '
  logM_TA_range           = '10.0,14.0,    '

  logA200_label           = '$\log_{10} A_{200}$,  '
  logA200_range           = '-2.0,2.0,     '

# Now parse target filename and select parameters to plot:  
  
  stem = string.split(filename,'_')[1]
  
  if stem == 'pairs':

    Npars = 8
    hline1 = '# importance, ' + logM_MW_label + logM_M31_label \
                     + D_M31_label + vr_M31_label + logMratio_label \
                     + logM_LG_pair_label \
                     + logM_TA_label + logA200_label 
                     # + vt_M31_label
    hline2 = '# 0,1,        ' + logM_MW_range + logM_M31_range \
                     + D_M31_range + vr_M31_range + logMratio_range \
                     + logM_LG_pair_range \
                     + logM_TA_range + logA200_range
                     # + vt_M31_range
    
    outbundle = numpy.zeros([n_given,Npars+1])
    outbundle[:,0] = w[index]
    outbundle[:,1] = numpy.log10(MW.M200[index])
    outbundle[:,2] = numpy.log10(M31.M200[index])
    outbundle[:,3] = M31.D[index]*1000.0
    outbundle[:,4] = M31.vr[index]
    outbundle[:,5] = numpy.log10(M31.M200[index]/MW.M200[index])
    outbundle[:,6] = numpy.log10(MW.M200[index]+M31.M200[index])
    outbundle[:,7] = numpy.log10(M31.M_TA[index])
    outbundle[:,8] = numpy.log10(M31.M_TA[index]/(MW.M200[index]+M31.M200[index]))
    # outbundle[:,9] = M31.vt[index]
    
  elif stem == 'triplets':

    Npars = 13
    hline1 = '# importance, ' \
                     + logM_MW_label + logM_M31_label + logM_M33_label \
                     + D_M31_label + vr_M31_label \
                     + D_M33_label + vr_M33_label \
                     + logMratio_label \
                     + logM_LG_pair_label + logM_LG_triplet_label \
                     + angsep_label \
                     + logM_TA_label + logA200_label
                     # + vt_M31_label
    hline2 = '# 0,1,        ' \
                     + logM_MW_range + logM_M31_range + logM_M33_range \
                     + D_M31_range + vr_M31_range \
                     + D_M33_range + vr_M33_range \
                     + logMratio_range \
                     + logM_LG_pair_range + logM_LG_triplet_range \
                     + angsep_range \
                     + logM_TA_range + logA200_range
                     # + vt_M31_range

    outbundle = numpy.zeros([n_given,Npars+1])
    outbundle[:,0] = w[index]
    outbundle[:,1] = numpy.log10(MW.M200[index])
    outbundle[:,2] = numpy.log10(M31.M200[index])
    outbundle[:,3] = numpy.log10(M33.M200[index])
    outbundle[:,4] = M31.D[index]*1000.0
    outbundle[:,5] = M31.vr[index]
    outbundle[:,6] = M33.D[index]*1000.0
    outbundle[:,7] = M33.vr[index]
    outbundle[:,8] = numpy.log10(M31.M200[index]/MW.M200[index])
    outbundle[:,9] = numpy.log10(MW.M200[index]+M31.M200[index])
    outbundle[:,10] = numpy.log10(MW.M200[index]+M31.M200[index]+M33.M200[index])
    outbundle[:,11] = M33.angsep[index]
    outbundle[:,12] = numpy.log10(M31.M_TA[index])
    outbundle[:,13] = numpy.log10(M31.M_TA[index]/(MW.M200[index]+M31.M200[index]))
    # outbundle[:,14] = M31.vt[index]
    
  else:
  
    print "Error: unrecognised target filename: "+filename

# # Now trim outbundle to random selection of n lines:
#   numpy.random.shuffle(outbundle)
#   outbundle = outbundle[0:n_required,:]

# # Now trim outbundle to (last) n halos with highest weight:
#   outbundle.sort(0)
#   outbundle = outbundle[n_given-n_required:n_given,:]

# The actual writing:

  # This is what I want to do - but header is not implemented in my 
  # numpy :-(  
  #   hdr = hline1+'\n'+hline2
  #   numpy.savetxt(filename, outbundle, header=hdr)

  output = open(filename,'w')
  output.write("%s\n" % hline1)
  output.write("%s\n" % hline2)
  output.close()
  numpy.savetxt('junk', outbundle)
  cat = subprocess.call("cat junk >> " + filename, shell=True)
  rm = subprocess.call("rm junk", shell=True)
  if cat != 0 or rm != 0:
    print "Error: write subprocesses failed in some way :-/"
    sys.exit()

  if vb: print "Written",n,"lines to",filename

  return
  
# ======================================================================

# NB. Files are named ".fits" but are actually of type tabletype, which 
# may well not be FITS...

def write_fits_files(kind,index,MW,M31,M33,w,vb):

  filename = "LG_"+kind+"_MW.fits"
  t = MW.rows(index[0])
  t.add_column('w',w)
  t.write(filename,type=tabletype,overwrite=True)
  if vb: print "Written",len(t),"rows to",filename
  
  filename = "LG_"+kind+"_M31.fits"
  t = M31.rows(index[0])
  t.add_column('w',w)
  t.write(filename,type=tabletype,overwrite=True)
  if vb: print "Written",len(t),"rows to",filename
  
  if kind == 'triplets':
     filename = "LG_"+kind+"_M33.fits"
     t = M33.rows(index[0])
     t.add_column('w',w)
     t.write(filename,type=tabletype,overwrite=True)
     if vb: print "Written",len(t),"rows to",filename
  
  return
  
# ======================================================================

# If called as a script, the python variable __name__ will be set to 
# "__main__" - so test for this, and execute the main program if so.
# Writing it like this allows the function plot_oguri_lens to be called
# from the python command line as well as from the unix prompt.

if __name__ == '__main__':
  LGHPinfer(sys.argv)

# ======================================================================


