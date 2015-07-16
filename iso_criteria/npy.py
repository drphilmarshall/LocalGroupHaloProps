import numpy as np
import os
import sys
sys.path.append('/u/ki/yymao/pyscripts')
from helpers.SimulationAnalysis import readHlist

#base_path = '/nfs/slac/g/ki/ki20/cosmo/behroozi/Consuelo/'
base_path = '/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/BGC2/'

box = int(sys.argv[1])

#os.mkdir('/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/'+str(box))
#halos = readHlist(base_path+str(box)+'/hlists/hlist_1.00000.list')
halos = readHlist(base_path+str(box)+'/'+str(box)+'bgc2_z000.txt')

np.save('/afs/slac.stanford.edu/u/ki/mwillia1/Thesis/LocalGroupHaloProps/iso_criteria/'+str(box)+'bgc2.npy', halos)
#np.save('/lustre/ki/pfs/mwillia1/LG_project/Consuelo_Boxes/'+str(box)+'/'+str(box)+'hlist', halos)
